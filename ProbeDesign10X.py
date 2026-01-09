"""
Probe design for 10x Genomics Visium and Chromium Flex custom probes.
Based on 10x Genomics technical note CG000621 Rev D:
https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import io
import json
import os
import re
import tempfile
import time
import urllib.parse
import urllib.request
import pandas as pd


# ---------------------------
# Helpers
# ---------------------------

_DNA_COMP = str.maketrans("ACGT", "TGCA")


def _clean_dna(seq: str) -> str:
    seq = re.sub(r"\s+", "", seq).upper().replace("U", "T")
    if not re.fullmatch(r"[ACGTN]+", seq):
        raise ValueError("Sequence contains non-ACGTN characters (U is allowed).")
    return seq


def revcomp(seq: str) -> str:
    seq = _clean_dna(seq)
    return seq.translate(_DNA_COMP)[::-1]


def read_fasta_records(path: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    seq_parts: list[str] = []
    with open(path, "r", encoding="ascii") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, _clean_dna("".join(seq_parts))))
                header = line[1:].strip() or f"seq_{len(records)+1:03d}"
                seq_parts = []
                continue
            seq_parts.append(line)
    if header is not None:
        records.append((header, _clean_dna("".join(seq_parts))))
    if not records:
        raise ValueError(f"No sequences found in FASTA: {path}")
    return records


def read_single_fasta(path: str) -> str:
    records = read_fasta_records(path)
    if len(records) != 1:
        raise ValueError(f"Expected 1 FASTA record, found {len(records)} in {path}.")
    return records[0][1]


def _gc_frac(s: str) -> float:
    gc = sum(1 for c in s if c in "GC")
    return gc / len(s) if s else 0.0


def _max_homopolymer_run(s: str) -> int:
    best = cur = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


# ---------------------------
# Chemistry constants (CG000621 Rev D)
# ---------------------------

VISIUM_LHS_PREFIX = "CCTTGGCACCCGAGAATTCCA"
VISIUM_RHS_POLYA = "A" * 30  # RHS has 30 A tail  [oai_citation:2‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)

FLEX_RHS_CONST5 = "ACGCGGTTAGCACGTA"
FLEX_RHS_CONST3 = "CGGTCCTAGCAA"
FLEX_SINGLEPLEX_BC001 = "ACTTTAGG"

FLEX_MULTIPLEX_BARCODES = {
    "BC001": "ACTTTAGG",
    "BC002": "AACGGGAA",
    "BC003": "AGTAGGCT",
    "BC004": "ATGTTGAC",
    "BC005": "ACAGACCT",
    "BC006": "ATCCCAAC",
    "BC007": "AAGTAGAG",
    "BC008": "AGCTGTGA",
    "BC009": "ACAGTCTG",
    "BC010": "AGTGAGTG",
    "BC011": "AGAGGCAA",
    "BC012": "ACTACTCA",
    "BC013": "ATACGTCA",
    "BC014": "ATCATGTG",
    "BC015": "AACGCCGA",
    "BC016": "ATTCGGTT",
}  #  [oai_citation:3‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)

FLEX_V2_RHS_PCONST = "CCCATATAAGAAA"
FLEX_V2_RHS_PCS1 = "CGGTCCTAGCAA"


# ---------------------------
# Probe design + formatting
# ---------------------------

@dataclass(frozen=True)
class Probe:
    start: int         # 0-based on input (sense / mRNA-like) sequence
    target50: str      # 50 nt window on input sequence (sense)
    seg1: str          # first 25 nt on input (sense)
    seg2: str          # last 25 nt on input (sense)
    lhs_target: str    # probe target_LHS portion (revcomp(seg1))
    rhs_target: str    # probe target_RHS portion (revcomp(seg2))
    gc_lhs: float
    gc_rhs: float
    max_hpoly: int
    score: float
    offtarget_hits: int = 0
    offtarget_subjects: tuple[str, ...] = ()


def _write_fasta(records: list[tuple[str, str]], path: str) -> None:
    with open(path, "w", encoding="ascii") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def _run_blast_remote_ncbi(
    fasta_text: str,
    blastdb: str,
    entrez_query: str | None,
    poll_interval: float,
    timeout_s: float,
    verbose: bool,
    cache: dict[str, str] | None,
    cache_key: str | None,
) -> pd.DataFrame:
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": blastdb,
        "QUERY": fasta_text,
        "FILTER": "F",
        "FORMAT_TYPE": "Tabular",
    }
    if entrez_query:
        params["ENTREZ_QUERY"] = entrez_query
    if cache is not None and cache_key and cache_key in cache:
        result_text = cache[cache_key]
        if verbose:
            print(f"[blast] cache hit for {blastdb}", flush=True)
        return _parse_blast_table(result_text)

    data = urllib.parse.urlencode(params).encode("ascii")
    if verbose:
        detail = f" db={blastdb}"
        if entrez_query:
            detail += f" entrez_query={entrez_query}"
        print(f"[offtarget] submit BLAST query:{detail}", flush=True)
    with urllib.request.urlopen("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=data) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    rid_match = re.search(r"RID = ([A-Z0-9-]+)", text)
    if not rid_match:
        raise RuntimeError(f"NCBI BLAST did not return RID. Response: {text[:500]}")
    rid = rid_match.group(1)

    start = time.time()
    if verbose:
        print(f"[blast] RID={rid} submitted", flush=True)
    while True:
        if time.time() - start > timeout_s:
            raise RuntimeError(f"NCBI BLAST timed out after {timeout_s} seconds for RID {rid}.")
        query = urllib.parse.urlencode({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"})
        with urllib.request.urlopen(f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?{query}") as resp:
            info = resp.read().decode("utf-8", errors="replace")
        if "Status=READY" in info:
            if "ThereAreHits=yes" not in info:
                if verbose:
                    print(f"[offtarget] RID={rid} complete, no hits", flush=True)
                return pd.DataFrame(
                    columns=[
                        "name",
                        "sseqid",
                        "pident",
                        "length",
                        "mismatch",
                        "gapopen",
                        "qstart",
                        "qend",
                        "sstart",
                        "send",
                        "evalue",
                        "bitscore",
                    ]
                )
            break
        if "Status=FAILED" in info:
            raise RuntimeError(f"NCBI BLAST failed for RID {rid}.")
        if verbose:
            print(f"[blast] RID={rid} still running...", flush=True)
        time.sleep(poll_interval)

    result_query = urllib.parse.urlencode(
        {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Tabular"}
    )
    if verbose:
        print(f"[blast] RID={rid} fetching results", flush=True)
    with urllib.request.urlopen(f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?{result_query}") as resp:
        result_text = resp.read().decode("utf-8", errors="replace")

    if cache is not None and cache_key:
        cache[cache_key] = result_text
    return _parse_blast_table(result_text)


def _parse_blast_table(result_text: str) -> pd.DataFrame:
    blast_columns = [
        "name",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    if not result_text.strip():
        return pd.DataFrame(columns=blast_columns)
    return pd.read_csv(io.StringIO(result_text), sep="\t", header=None, names=blast_columns)


def _load_blast_cache(path: str | None) -> dict[str, str]:
    if not path or not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _save_blast_cache(path: str | None, cache: dict[str, str]) -> None:
    if not path:
        return
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cache, fh, indent=2, sort_keys=True)


def _summarize_offtargets(
    blast_res: pd.DataFrame,
    subject_id_sep: str | None,
    target_id: str | None,
    min_mismatches: int,
) -> dict[str, tuple[int, tuple[str, ...]]]:
    if blast_res.empty:
        return {}

    res = blast_res.copy()
    if subject_id_sep:
        res["subject_id"] = res["sseqid"].str.split(subject_id_sep).str[0]
    else:
        res["subject_id"] = res["sseqid"]

    hits = res[res["mismatch"] < min_mismatches]
    if target_id is not None:
        hits = hits[hits["subject_id"] != target_id]

    summary = {}
    for name, subdf in hits.groupby("name"):
        subjects = tuple(sorted(set(subdf["subject_id"].tolist())))
        summary[name] = (len(subjects), subjects)
    return summary


def design_custom_probes(
    seq: str,
    assay: str,                         # "visium" or "flex"
    version: str | None = None,         # e.g. "v1"/"v2"/"hd" for visium; ignored for flex in this doc
    flex_mode: str | None = None,       # "singleplex" or "multiplex" (required if assay="flex")
    flex_version: str | None = None,    # "v1", "v2", or "v2_4plex" (only if assay="flex")
    flex_barcode: str = "BC001",        # for multiplex: BC001..BC016; for singleplex ignored
    probe_len: int = 50,
    split: int = 25,
    gc_min: float = 0.44,
    gc_max: float = 0.72,
    max_homopolymer: int = 6,
    min_spacing: int = 50,
    max_probes: int = 12,
    nn: str = "NN",                     # Flex RHS includes NN  [oai_citation:4‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)
    blast_for_offtarget: bool = False,
    blast_dbs: tuple[str, ...] = ("nt", "refseq_rna"),
    blast_organism: str | None = None,
    blast_min_mismatches: int = 5,
    blast_max_hits: int = 1,
    blast_subject_id_sep: str | None = None,
    blast_target_id: str | None = None,
    blast_poll_interval: float = 5.0,
    blast_timeout_s: float = 600.0,
    blast_verbose: bool = True,
    probe_name_prefix: str | None = None,
    ligation_requires_a: bool = True,
    blast_cache_path: str | None = None,
    auto_pick: bool = True,
) -> pd.DataFrame:
    """
    Designs split probes from an input (mRNA-like, 5'->3') sequence.
    Outputs order-ready oligos with the correct constant arms for Visium or Chromium Flex (CG000621 Rev D).

    Key design constraint for Visium: LHS probe's 3' base must be T, meaning seg1[-1] on the input must be A.  [oai_citation:5‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)

    Flex v2 (GEM-X) uses different RHS constant sequences and omits barcodes/NNs.
    Use flex_version="v2_4plex" to apply the pCS1 RHS for the 4-sample kit.

    Off-target screening (optional) uses NCBI remote BLAST. Set blast_for_offtarget=True
    and provide blast_organism (e.g., "human" or "mouse") to constrain results.
    If blast_target_id is None, blast_max_hits should usually be 1 to allow the primary hit.
    Use blast_cache_path to cache results by query + database.
    """
    assay = assay.lower()
    if assay not in {"visium", "flex"}:
        raise ValueError("assay must be 'visium' or 'flex'")

    if split != 25 or probe_len != 50:
        raise ValueError("This implementation assumes 50 nt probes split 25/25 (as in CG000621 guidance).")

    if assay == "flex":
        if flex_version is None:
            flex_version = "v1"
        flex_version = flex_version.lower()
        if flex_version not in {"v1", "v2", "v2_4plex"}:
            raise ValueError("For assay='flex', flex_version should be 'v1', 'v2', or 'v2_4plex'.")

        if flex_mode is None:
            flex_mode = "singleplex" if flex_version in {"v2", "v2_4plex"} else None
        if flex_mode is None or flex_mode.lower() not in {"singleplex", "multiplex"}:
            raise ValueError("For assay='flex', set flex_mode to 'singleplex' or 'multiplex'.")
        flex_mode = flex_mode.lower()

        if flex_version == "v1":
            if flex_mode == "multiplex" and flex_barcode not in FLEX_MULTIPLEX_BARCODES:
                raise ValueError("flex_barcode must be one of BC001..BC016 for multiplex.")
    else:
        # Visium v1/v2/HD share the same ordering structure in this tech note  [oai_citation:6‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)
        if version is not None:
            version = version.lower()
            if version not in {"v1", "v2", "hd"}:
                raise ValueError("For assay='visium', version should be one of: v1, v2, hd (or None).")

    seq = _clean_dna(seq)
    if len(seq) < probe_len:
        raise ValueError(f"Sequence length {len(seq)} is shorter than probe_len {probe_len}.")

    if auto_pick:
        max_probes = 3
        min_spacing = max(min_spacing, probe_len)

    cands: list[Probe] = []
    for start in range(0, len(seq) - probe_len + 1):
        window = seq[start : start + probe_len]
        if "N" in window:
            continue

        seg1 = window[:split]
        seg2 = window[split:]

        lhs_target = revcomp(seg1)
        rhs_target = revcomp(seg2)

        # Ligation-junction constraint: LHS 3' base (25th in probe) must be T.  [oai_citation:7‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)
        if ligation_requires_a and lhs_target[-1] != "T":
            continue

        gc_l = _gc_frac(lhs_target)
        gc_r = _gc_frac(rhs_target)
        if not (gc_min <= gc_l <= gc_max and gc_min <= gc_r <= gc_max):
            continue

        hp = _max_homopolymer_run(lhs_target + rhs_target)
        if hp > max_homopolymer:
            continue

        if auto_pick:
            target_gc = 0.5
            gc_weight = 2.0
        else:
            target_gc = (gc_min + gc_max) / 2
            gc_weight = 1.0
        gc_pen = gc_weight * (abs(gc_l - target_gc) + abs(gc_r - target_gc))
        score = 1.0 - gc_pen - 0.05 * max(0, hp - 3)

        cands.append(
            Probe(
                start=start,
                target50=window,
                seg1=seg1,
                seg2=seg2,
                lhs_target=lhs_target,
                rhs_target=rhs_target,
                gc_lhs=gc_l,
                gc_rhs=gc_r,
                max_hpoly=hp,
                score=score,
            )
        )

    organism_query = None
    if blast_organism:
        org_map = {
            "human": "Homo sapiens[Organism]",
            "mouse": "Mus musculus[Organism]",
        }
        organism_query = org_map.get(blast_organism.lower(), blast_organism)

    if blast_for_offtarget:
        cache = _load_blast_cache(blast_cache_path)
        cache_dirty = False
        records = [(f"cand_{i}", p.target50) for i, p in enumerate(cands)]
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "queries.fa")
            _write_fasta(records, fasta_path)
            fasta_text = open(fasta_path, "r", encoding="ascii").read()
            blast_frames = []
            if blast_for_offtarget:
                for db in blast_dbs:
                    cache_key = None
                    if blast_cache_path:
                        key_base = f"{db}|{organism_query or ''}|{fasta_text}"
                        cache_key = hashlib.sha256(key_base.encode("utf-8")).hexdigest()
                        had_in_cache = cache_key in cache
                    blast_frames.append(
                        _run_blast_remote_ncbi(
                            fasta_text=fasta_text,
                            blastdb=db,
                            entrez_query=organism_query,
                            poll_interval=blast_poll_interval,
                            timeout_s=blast_timeout_s,
                            verbose=blast_verbose,
                            cache=cache,
                            cache_key=cache_key,
                        )
                    )
                    if cache_key and not had_in_cache and cache_key in cache:
                        cache_dirty = True
            blast_res = pd.concat(blast_frames, ignore_index=True)
        if blast_cache_path and cache_dirty:
            _save_blast_cache(blast_cache_path, cache)
        summary = _summarize_offtargets(
            blast_res,
            subject_id_sep=blast_subject_id_sep,
            target_id=blast_target_id,
            min_mismatches=blast_min_mismatches,
        )
        updated = []
        for i, p in enumerate(cands):
            name = f"cand_{i}"
            hits, subjects = summary.get(name, (0, ()))
            updated.append(
                replace(p, offtarget_hits=hits, offtarget_subjects=subjects)
            )
        cands = [p for p in updated if p.offtarget_hits <= blast_max_hits]

    # greedy non-overlapping selection
    cands.sort(key=lambda p: p.score, reverse=True)
    selected: list[Probe] = []
    for p in cands:
        if len(selected) >= max_probes:
            break
        if all(abs(p.start - q.start) >= min_spacing for q in selected):
            selected.append(p)
    selected.sort(key=lambda p: p.start)

    # Format oligos for ordering per assay
    rows = []
    for i, p in enumerate(selected, start=1):
        if probe_name_prefix:
            probe_id = f"{probe_name_prefix}_{i:02d}"
        else:
            probe_id = f"probe_{i:02d}"

        if assay == "visium":
            # Visium ordering format  [oai_citation:8‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)
            lhs_oligo = f"{VISIUM_LHS_PREFIX}{p.lhs_target}"
            rhs_oligo = f"/5Phos/{p.rhs_target}{VISIUM_RHS_POLYA}"
            rhs_5phos = True
            rhs_note = "/5Phos/ required"
        else:
            lhs_oligo = f"{VISIUM_LHS_PREFIX}{p.lhs_target}"
            rhs_5phos = True
            if flex_version in {"v2", "v2_4plex"}:
                rhs_const = FLEX_V2_RHS_PCS1 if flex_version == "v2_4plex" else FLEX_V2_RHS_PCONST
                rhs_oligo = f"{p.rhs_target}{rhs_const}"
                rhs_note = "/5Phos/ required; GEM-X Flex v2"
                if flex_version == "v2_4plex":
                    rhs_note = "/5Phos/ required; GEM-X Flex v2 4-sample kit"
            else:
                # Flex ordering format (singleplex/multiplex)  [oai_citation:9‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)
                if flex_mode == "singleplex":
                    bc = FLEX_SINGLEPLEX_BC001
                else:
                    bc = FLEX_MULTIPLEX_BARCODES[flex_barcode]

                # target_RHS - CONST5 - NN - BARCODE - CONST3  [oai_citation:10‡Cloudinary](https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf)
                rhs_oligo = f"{p.rhs_target}{FLEX_RHS_CONST5}{nn}{bc}{FLEX_RHS_CONST3}"
                rhs_note = "/5Phos/ required; includes NN + probe barcode"
            if rhs_5phos:
                rhs_oligo = "/5Phos/" + rhs_oligo

        if assay == "visium":
            assay_version = version or "default"
        else:
            assay_version = flex_version or "v1"

        row = {
            "probe_id": probe_id,
            "assay": assay,
            "assay_version": assay_version,
            "start_0based": p.start,
            "target_window_50_sense": p.target50,
            "target_LHS_25_probe": p.lhs_target,
            "target_RHS_25_probe": p.rhs_target,
            "gc_LHS": p.gc_lhs,
            "gc_RHS": p.gc_rhs,
            "max_homopolymer": p.max_hpoly,
            "score": p.score,
            "LHS_order_5to3": lhs_oligo,
            "RHS_order_5to3": rhs_oligo,
        }
        if assay == "flex":
            row["flex_mode"] = flex_mode
            if flex_version == "v1" and flex_mode == "multiplex":
                row["flex_barcode"] = flex_barcode
        if blast_for_offtarget:
            row["offtarget_hits"] = p.offtarget_hits
            row["offtarget_subjects"] = ",".join(p.offtarget_subjects)

        rows.append(row)

    return pd.DataFrame(rows)


def design_custom_probes_from_fasta(
    fasta_path: str,
    *,
    assay: str,
    version: str | None = None,
    flex_mode: str | None = None,
    flex_version: str | None = None,
    flex_barcode: str = "BC001",
    probe_len: int = 50,
    split: int = 25,
    gc_min: float = 0.44,
    gc_max: float = 0.72,
    max_homopolymer: int = 6,
    min_spacing: int = 50,
    max_probes: int = 12,
    auto_pick: bool = True,
    nn: str = "NN",
    blast_for_offtarget: bool = False,
    blast_dbs: tuple[str, ...] = ("nt", "refseq_rna"),
    blast_organism: str | None = None,
    blast_min_mismatches: int = 5,
    blast_max_hits: int = 1,
    blast_subject_id_sep: str | None = None,
    blast_target_id: str | None = None,
    blast_poll_interval: float = 5.0,
    blast_timeout_s: float = 600.0,
    blast_verbose: bool = True,
    probe_name_prefix: str | None = None,
    ligation_requires_a: bool = True,
    blast_cache_path: str | None = None,
) -> pd.DataFrame:
    records = read_fasta_records(fasta_path)
    all_frames = []
    for name, seq in records:
        prefix = probe_name_prefix or name
        df = design_custom_probes(
            seq=seq,
            assay=assay,
            version=version,
            flex_version=flex_version,
            flex_mode=flex_mode,
            flex_barcode=flex_barcode,
            probe_len=probe_len,
            split=split,
            gc_min=gc_min,
            gc_max=gc_max,
            max_homopolymer=max_homopolymer,
            min_spacing=min_spacing,
            max_probes=max_probes,
            auto_pick=auto_pick,
            nn=nn,
            blast_for_offtarget=blast_for_offtarget,
            blast_dbs=blast_dbs,
            blast_organism=blast_organism,
            blast_min_mismatches=blast_min_mismatches,
            blast_max_hits=blast_max_hits,
            blast_subject_id_sep=blast_subject_id_sep,
            blast_target_id=blast_target_id,
            blast_poll_interval=blast_poll_interval,
            blast_timeout_s=blast_timeout_s,
            blast_verbose=blast_verbose,
            probe_name_prefix=prefix,
            ligation_requires_a=ligation_requires_a,
            blast_cache_path=blast_cache_path,
        )
        if len(records) > 1:
            df.insert(0, "sequence_id", name)
        all_frames.append(df)

    return pd.concat(all_frames, ignore_index=True) if all_frames else pd.DataFrame()


def plot_probe_hybridization(
    probes_df: pd.DataFrame,
    seq_len: int | None = None,
    max_probes_to_plot: int = 100,
    sort_by: str = "start_0based",
    figsize: tuple[float, float] = (10.0, 4.0),
):
    """
    Plot probe hybridization positions along a reference sequence using Matplotlib.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("Matplotlib is required for plotting.") from exc

    if probes_df.empty:
        raise ValueError("probes_df is empty.")
    if sort_by not in probes_df.columns:
        raise ValueError(f"sort_by column '{sort_by}' not found in probes_df.")
    if "start_0based" not in probes_df.columns:
        raise ValueError("probes_df must include 'start_0based'.")

    df = probes_df.sort_values(sort_by).reset_index(drop=True)
    if len(df) > max_probes_to_plot:
        df = df.head(max_probes_to_plot)

    if seq_len is None:
        if "target_window_50_sense" in df.columns:
            win_len = df["target_window_50_sense"].str.len().max()
        else:
            win_len = 50
        seq_len = int(df["start_0based"].max() + win_len)

    fig, ax = plt.subplots(figsize=figsize)
    y_positions = range(len(df))
    starts = df["start_0based"].astype(int).tolist()
    if "target_window_50_sense" in df.columns:
        lengths = df["target_window_50_sense"].str.len().fillna(50).astype(int).tolist()
    else:
        lengths = [50] * len(df)

    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#1f77b4"])
    colors = [cycle[i % len(cycle)] for i in range(len(df))]
    ax.barh(list(y_positions), lengths, left=starts, height=0.6, color=colors)
    ax.set_xlim(0, seq_len)
    ax.set_xlabel("Position on reference (0-based)")
    ax.set_ylabel("Probe index")
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(df["probe_id"].tolist() if "probe_id" in df.columns else [str(i) for i in y_positions])
    ax.set_title("Probe hybridization layout")
    ax.invert_yaxis()
    fig.tight_layout()
    return ax


def to_idt_opools(
    probes_df: pd.DataFrame,
    pool_name: str = "poolOne",
) -> pd.DataFrame:
    """
    Convert probe output into IDT oPools order format.
    """
    if probes_df.empty:
        raise ValueError("probes_df is empty.")
    if "LHS_order_5to3" not in probes_df.columns or "RHS_order_5to3" not in probes_df.columns:
        raise ValueError("probes_df must include LHS_order_5to3 and RHS_order_5to3.")

    sequences = []
    sequences.extend(probes_df["LHS_order_5to3"].tolist())
    sequences.extend(probes_df["RHS_order_5to3"].tolist())
    return pd.DataFrame({"Pool name": [pool_name] * len(sequences), "Sequence": sequences})


# ---------------------------
# Example
# ---------------------------
# seq = read_single_fasta("my_construct.fa")
# visium_df = design_custom_probes(seq, assay="visium", version="hd", max_probes=12)
# visium_df = design_custom_probes(
#     seq,
#     assay="visium",
#     version="hd",
#     blast_for_offtarget=True,
#     blast_organism="human",
#     blast_target_id="MY_GENE",
#     blast_subject_id_sep="::",
#     blast_max_hits=0,
#     blast_cache_path="blast_cache.json",
# )
# flex_df   = design_custom_probes(seq, assay="flex", flex_mode="multiplex", flex_barcode="BC003", max_probes=12)
# flex_v2_df = design_custom_probes(
#     seq,
#     assay="flex",
#     flex_version="v2",
#     flex_mode="singleplex",
#     max_probes=12,
# )
# flex_v2_4plex_df = design_custom_probes(
#     seq,
#     assay="flex",
#     flex_version="v2_4plex",
#     flex_mode="singleplex",
#     max_probes=12,
# )
# visium_df.to_csv("custom_probes_visium.csv", index=False)
# flex_df.to_csv("custom_probes_flex.csv", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Design 10x custom probes.")
    parser.add_argument("fasta", help="Input FASTA with a single sequence.")
    parser.add_argument("--assay", required=True, choices=["visium", "flex"])
    parser.add_argument("--version", choices=["v1", "v2", "hd"])
    parser.add_argument("--flex-version", choices=["v1", "v2", "v2_4plex"])
    parser.add_argument("--flex-mode", choices=["singleplex", "multiplex"])
    parser.add_argument("--flex-barcode", default="BC001")
    parser.add_argument("--probe-name-prefix")
    parser.add_argument("--max-probes", type=int, default=12)
    parser.add_argument("--min-spacing", type=int, default=50)
    parser.add_argument("--no-ligation-constraint", action="store_true")
    parser.add_argument("--blast", action="store_true")
    parser.add_argument("--blast-dbs", default="nt,refseq_rna")
    parser.add_argument("--blast-organism")
    parser.add_argument("--blast-min-mismatches", type=int, default=5)
    parser.add_argument("--blast-max-hits", type=int, default=1)
    parser.add_argument("--blast-subject-id-sep")
    parser.add_argument("--blast-target-id")
    parser.add_argument("--blast-poll-interval", type=float, default=5.0)
    parser.add_argument("--blast-timeout-s", type=float, default=600.0)
    parser.add_argument("--blast-quiet", action="store_true")
    parser.add_argument("--blast-cache")
    parser.add_argument("--out")
    parser.add_argument("--plot")

    args = parser.parse_args()
    out_df = design_custom_probes_from_fasta(
        args.fasta,
        assay=args.assay,
        version=args.version,
        flex_version=args.flex_version,
        flex_mode=args.flex_mode,
        flex_barcode=args.flex_barcode,
        max_probes=args.max_probes,
        min_spacing=args.min_spacing,
        blast_for_offtarget=args.blast,
        blast_dbs=tuple(x for x in args.blast_dbs.split(",") if x),
        blast_organism=args.blast_organism,
        blast_min_mismatches=args.blast_min_mismatches,
        blast_max_hits=args.blast_max_hits,
        blast_subject_id_sep=args.blast_subject_id_sep,
        blast_target_id=args.blast_target_id,
        blast_poll_interval=args.blast_poll_interval,
        blast_timeout_s=args.blast_timeout_s,
        blast_verbose=not args.blast_quiet,
        probe_name_prefix=args.probe_name_prefix,
        ligation_requires_a=not args.no_ligation_constraint,
        blast_cache_path=args.blast_cache,
    )

    if args.out:
        out_df.to_csv(args.out, index=False)
    else:
        print(out_df.to_csv(index=False).strip())

    if args.plot:
        records = read_fasta_records(args.fasta)
        if len(records) != 1:
            raise ValueError("--plot requires a single FASTA record.")
        ax = plot_probe_hybridization(out_df, seq_len=len(records[0][1]))
        ax.figure.savefig(args.plot, dpi=200)
