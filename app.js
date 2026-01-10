// -------- core utilities --------

function parseFasta(text) {
  const records = [];
  let header = null;
  let seq = [];
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith(">")) {
      if (header) records.push({ id: header, seq: seq.join("") });
      header = line.slice(1).trim().split(/\s+/)[0] || `seq${records.length + 1}`;
      seq = [];
    } else {
      seq.push(line.replace(/\s+/g, ""));
    }
  }
  if (header) records.push({ id: header, seq: seq.join("") });
  return records;
}

function cleanDna(seq) {
  const s = seq.replace(/\s+/g, "").toUpperCase().replace(/U/g, "T");
  if (!s) return "";
  if (/[^ACGTN]/.test(s)) {
    throw new Error("Sequence contains non-ACGTN characters (U allowed).");
  }
  return s;
}

const DNA_COMP = {
  A: "T",
  C: "G",
  G: "C",
  T: "A",
  N: "N",
};

function revcomp(seq) {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) {
    out += DNA_COMP[seq[i]] || "N";
  }
  return out;
}

function gcFrac(s) {
  let gc = 0;
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if (c === "G" || c === "C") gc++;
  }
  return gc / s.length;
}

function maxHomopolymerRun(s) {
  let best = 1, cur = 1;
  for (let i = 1; i < s.length; i++) {
    if (s[i] === s[i - 1]) {
      cur++;
      if (cur > best) best = cur;
    } else {
      cur = 1;
    }
  }
  return best;
}

function padNumber(n, width) {
  const s = String(n);
  return s.length >= width ? s : "0".repeat(width - s.length) + s;
}

function downloadText(filename, text) {
  const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function toCsv(rows, columns) {
  const esc = (v) => {
    const s = v == null ? "" : String(v);
    if (/[",\n]/.test(s)) return `"${s.replace(/"/g, '""')}"`;
    return s;
  };
  const header = columns.map(esc).join(",");
  const lines = rows.map(r => columns.map(c => esc(r[c])).join(","));
  return [header, ...lines].join("\n");
}

// -------- design logic --------

const PROBE_LEN = 50;
const SPLIT_LEN = 25;

const VISIUM_LHS_PREFIX = "CCTTGGCACCCGAGAATTCCA";
const VISIUM_RHS_POLYA = "A".repeat(30);

const FLEX_RHS_CONST5 = "ACGCGGTTAGCACGTA";
const FLEX_RHS_CONST3 = "CGGTCCTAGCAA";
const FLEX_SINGLEPLEX_BC001 = "ACTTTAGG";
const FLEX_MULTIPLEX_BARCODES = {
  BC001: "ACTTTAGG",
  BC002: "AACGGGAA",
  BC003: "AGTAGGCT",
  BC004: "ATGTTGAC",
  BC005: "ACAGACCT",
  BC006: "ATCCCAAC",
  BC007: "AAGTAGAG",
  BC008: "AGCTGTGA",
  BC009: "ACAGTCTG",
  BC010: "AGTGAGTG",
  BC011: "AGAGGCAA",
  BC012: "ACTACTCA",
  BC013: "ATACGTCA",
  BC014: "ATCATGTG",
  BC015: "AACGCCGA",
  BC016: "ATTCGGTT",
};

const FLEX_V2_RHS_PCONST = "CCCATATAAGAAA";
const FLEX_V2_RHS_PCS1 = "CGGTCCTAGCAA";

function buildCandidates(seq, params) {
  const candidates = [];
  for (let start = 0; start + PROBE_LEN <= seq.length; start++) {
    const window = seq.slice(start, start + PROBE_LEN);
    if (window.includes("N")) continue;

    const seg1 = window.slice(0, SPLIT_LEN);
    const seg2 = window.slice(SPLIT_LEN);

    const lhsTarget = revcomp(seg1);
    const rhsTarget = revcomp(seg2);

    if (params.ligationRequiresT && lhsTarget[lhsTarget.length - 1] !== "T") continue;

    const gcL = gcFrac(lhsTarget);
    const gcR = gcFrac(rhsTarget);
    if (gcL < params.gcMin || gcL > params.gcMax) continue;
    if (gcR < params.gcMin || gcR > params.gcMax) continue;

    const hpoly = maxHomopolymerRun(lhsTarget + rhsTarget);
    if (hpoly > params.maxHpoly) continue;

    const hasAutoPick = Number.isFinite(params.autoPickCount);
    const targetGc = hasAutoPick ? 0.5 : (params.gcMin + params.gcMax) / 2;
    const gcWeight = hasAutoPick ? 2.0 : 1.0;
    const gcPen = gcWeight * (Math.abs(gcL - targetGc) + Math.abs(gcR - targetGc));
    const score = 1.0 - gcPen - 0.05 * Math.max(0, hpoly - 3);

    candidates.push({
      start,
      target50: window,
      seg1,
      seg2,
      lhsTarget,
      rhsTarget,
      gcL,
      gcR,
      maxHpoly: hpoly,
      score,
    });
  }

  return candidates;
}

function selectWithSpacing(cands, minSpacing, maxProbes) {
  const spacing = Math.max(0, minSpacing);
  const limit = maxProbes > 0 ? maxProbes : Infinity;
  const sorted = cands.slice().sort((a, b) => b.score - a.score);
  const selected = [];
  for (const cand of sorted) {
    if (selected.length >= limit) break;
    const ok = selected.every(ch => Math.abs(cand.start - ch.start) >= spacing);
    if (!ok) continue;
    selected.push(cand);
  }
  return selected.sort((a, b) => a.start - b.start);
}

function formatProbeRow(recId, probeIndex, cand, params) {
  const probeId = `${recId}_${padNumber(probeIndex, 2)}`;

  const row = {
    sequence_id: recId,
    probe_id: probeId,
    assay: params.assay,
    assay_version: params.assay === "visium" ? (params.visiumVersion || "default") : params.flexVersion,
    start_0based: cand.start,
    target_window_50_sense: cand.target50,
    target_LHS_25_probe: cand.lhsTarget,
    target_RHS_25_probe: cand.rhsTarget,
    gc_LHS: +cand.gcL.toFixed(4),
    gc_RHS: +cand.gcR.toFixed(4),
    max_homopolymer: cand.maxHpoly,
    score: +cand.score.toFixed(4),
    LHS_order_5to3: "",
    RHS_order_5to3: "",
    flex_mode: "",
    flex_barcode: "",
  };

  row.LHS_order_5to3 = `${VISIUM_LHS_PREFIX}${cand.lhsTarget}`;

  if (params.assay === "visium") {
    row.RHS_order_5to3 = `/5Phos/${cand.rhsTarget}${VISIUM_RHS_POLYA}`;
  } else {
    row.flex_mode = params.flexMode;
    if (params.flexVersion === "v2" || params.flexVersion === "v2_4plex") {
      const rhsConst = params.flexVersion === "v2_4plex" ? FLEX_V2_RHS_PCS1 : FLEX_V2_RHS_PCONST;
      row.RHS_order_5to3 = `/5Phos/${cand.rhsTarget}${rhsConst}`;
    } else {
      const barcode = params.flexMode === "multiplex"
        ? FLEX_MULTIPLEX_BARCODES[params.flexBarcode]
        : FLEX_SINGLEPLEX_BC001;
      row.flex_barcode = params.flexMode === "multiplex" ? params.flexBarcode : "";
      row.RHS_order_5to3 = `/5Phos/${cand.rhsTarget}${FLEX_RHS_CONST5}${params.nn}${barcode}${FLEX_RHS_CONST3}`;
    }
  }

  return row;
}

function designProbes(fastaText, params) {
  const records = parseFasta(fastaText);
  if (!records.length) {
    throw new Error("No FASTA records found.");
  }

  const stats = {
    targets: records.length,
    total_bases: 0,
    candidates: 0,
    selected: 0,
  };

  const rows = [];
  for (const rec of records) {
    const seq = cleanDna(rec.seq);
    if (!seq) continue;
    if (seq.length < PROBE_LEN) {
      throw new Error(`Sequence ${rec.id} is shorter than ${PROBE_LEN} nt.`);
    }
    stats.total_bases += seq.length;

    const hasAutoPick = Number.isFinite(params.autoPickCount);
    const effectiveMinSpacing = hasAutoPick ? Math.max(params.minSpacing, PROBE_LEN) : params.minSpacing;
    const effectiveMaxProbes = hasAutoPick ? params.autoPickCount : Infinity;

    const candidates = buildCandidates(seq, params);
    stats.candidates += candidates.length;
    const selected = selectWithSpacing(candidates, effectiveMinSpacing, effectiveMaxProbes);
    stats.selected += selected.length;

    selected.forEach((cand, idx) => {
      rows.push(formatProbeRow(rec.id, idx + 1, cand, params));
    });
  }

  return { rows, stats };
}

// -------- table rendering --------

function renderTable(tableEl, rows, columns, sortState) {
  const thead = tableEl.querySelector("thead");
  const tbody = tableEl.querySelector("tbody");

  // header
  thead.innerHTML = "";
  const trh = document.createElement("tr");
  for (const col of columns) {
    const th = document.createElement("th");
    th.textContent = col + (sortState.col === col ? (sortState.dir === "asc" ? " ▲" : " ▼") : "");
    th.dataset.col = col;
    trh.appendChild(th);
  }
  thead.appendChild(trh);

  // body
  tbody.innerHTML = "";
  const frag = document.createDocumentFragment();
  for (const r of rows) {
    const tr = document.createElement("tr");
    for (const col of columns) {
      const td = document.createElement("td");
      const v = r[col];
      td.textContent = v == null ? "" : String(v);
      if (col.includes("order") || col.includes("target_") || col.includes("window_")) {
        td.classList.add("seq");
      }
      tr.appendChild(td);
    }
    frag.appendChild(tr);
  }
  tbody.appendChild(frag);
}

function sortRows(rows, col, dir) {
  const copy = rows.slice();
  const mul = dir === "asc" ? 1 : -1;
  copy.sort((a, b) => {
    const av = a[col], bv = b[col];
    const an = typeof av === "number" ? av : Number(av);
    const bn = typeof bv === "number" ? bv : Number(bv);
    const bothNum = Number.isFinite(an) && Number.isFinite(bn);

    if (bothNum) return (an - bn) * mul;
    return String(av ?? "").localeCompare(String(bv ?? "")) * mul;
  });
  return copy;
}

function filterRows(rows, query) {
  const q = query.trim().toLowerCase();
  if (!q) return rows;
  return rows.filter(r => {
    for (const v of Object.values(r)) {
      if (v == null) continue;
      if (String(v).toLowerCase().includes(q)) return true;
    }
    return false;
  });
}

// -------- oPools export --------

function toOpoolsRows(designRows, poolName) {
  // IDT oPools can be as simple as Name,Sequence (commonly accepted).
  // If you want other fields (Scale, Purification, Plate, Well, etc.), add columns here.
  const rows = [];
  const name = poolName || "poolOne";
  for (const r of designRows) {
    rows.push({ "Pool name": name, Sequence: r.LHS_order_5to3 });
  }
  for (const r of designRows) {
    rows.push({ "Pool name": name, Sequence: r.RHS_order_5to3 });
  }
  return rows;
}

// -------- UI wiring --------

const els = {
  fasta: document.getElementById("fasta"),
  file: document.getElementById("file"),
  example: document.getElementById("example"),
  run: document.getElementById("run"),
  status: document.getElementById("status"),
  stats: document.getElementById("stats"),
  search: document.getElementById("search"),
  dlCsv: document.getElementById("dlCsv"),
  dlOpools: document.getElementById("dlOpools"),
  table: document.getElementById("table"),
  poolName: document.getElementById("poolName"),
  assay: document.getElementById("assay"),
  visiumOptions: document.getElementById("visiumOptions"),
  visiumVersion: document.getElementById("visiumVersion"),
  flexOptions: document.getElementById("flexOptions"),
  flexVersion: document.getElementById("flexVersion"),
  flexMode: document.getElementById("flexMode"),
  flexBarcode: document.getElementById("flexBarcode"),
  flexNN: document.getElementById("flexNN"),
  probeLen: document.getElementById("probeLen"),
  minSpacing: document.getElementById("minSpacing"),
  gcMin: document.getElementById("gcMin"),
  gcMax: document.getElementById("gcMax"),
  maxHpoly: document.getElementById("maxHpoly"),
  autoPickCount: document.getElementById("autoPickCount"),
  ligationA: document.getElementById("ligationA"),
};

let lastDesign = { rows: [], stats: null };
let currentViewRows = [];
let sortState = { col: "probe_id", dir: "asc" };

// Columns for the “pandas-like” table
const TABLE_COLUMNS = [
  "sequence_id",
  "probe_id",
  "assay",
  "assay_version",
  "start_0based",
  "target_window_50_sense",
  "target_LHS_25_probe",
  "target_RHS_25_probe",
  "gc_LHS",
  "gc_RHS",
  "max_homopolymer",
  "score",
  "LHS_order_5to3",
  "RHS_order_5to3",
  "flex_mode",
  "flex_barcode",
];

function getParams() {
  const assay = els.assay.value;
  const visiumVersion = els.visiumVersion.value || "";
  const flexVersion = els.flexVersion.value;
  const flexMode = els.flexMode.value;
  const flexBarcode = els.flexBarcode.value;
  const nn = els.flexNN.value.trim().toUpperCase();

  if (assay === "visium" && visiumVersion && visiumVersion !== "v1/v2/HD") {
    throw new Error("Visium version must be v1/v2/HD.");
  }

  if (assay === "flex") {
    if (!["v1", "v2", "v2_4plex"].includes(flexVersion)) {
      throw new Error("Flex version must be v1, v2, or v2_4plex.");
    }
    if (flexVersion === "v1" && !["singleplex", "multiplex"].includes(flexMode)) {
      throw new Error("Flex mode must be singleplex or multiplex.");
    }
    if (flexVersion === "v1" && flexMode === "multiplex" && !(flexBarcode in FLEX_MULTIPLEX_BARCODES)) {
      throw new Error("Flex barcode must be BC001..BC016 for multiplex.");
    }
    if (flexVersion === "v1") {
      if (!/^[ACGTN]{2}$/.test(nn)) {
        throw new Error("Flex NN must be exactly 2 bases (A/C/G/T/N).");
      }
    }
  }

  const gcMin = Number(els.gcMin.value);
  const gcMax = Number(els.gcMax.value);
  if (!Number.isFinite(gcMin) || !Number.isFinite(gcMax) || gcMin < 0 || gcMax > 1 || gcMin > gcMax) {
    throw new Error("GC bounds must be numbers in [0,1] with GC min <= GC max.");
  }

  const minSpacing = Number(els.minSpacing.value);
  if (!Number.isFinite(minSpacing) || minSpacing < 0) {
    throw new Error("Min spacing must be a non-negative number.");
  }

  const autoPickRaw = els.autoPickCount.value.trim();
  let autoPickCount = null;
  if (autoPickRaw) {
    autoPickCount = Number(autoPickRaw);
    if (!Number.isFinite(autoPickCount) || autoPickCount < 1) {
      throw new Error("Auto-pick count must be >= 1, or blank for all.");
    }
  }

  const maxHpoly = Number(els.maxHpoly.value);
  if (!Number.isFinite(maxHpoly) || maxHpoly < 1) {
    throw new Error("Max homopolymer run must be >= 1.");
  }

  const probeLen = Number(els.probeLen.value);
  if (probeLen !== PROBE_LEN) {
    throw new Error(`Probe length must be ${PROBE_LEN} nt for 25/25 split.`);
  }

  return {
    assay,
    visiumVersion,
    flexVersion,
    flexMode,
    flexBarcode,
    nn: nn || "NN",
    probeLen,
    minSpacing,
    gcMin,
    gcMax,
    maxHpoly,
    autoPickCount,
    ligationRequiresT: els.ligationA.checked,
  };
}

function updateAssayControls() {
  const isVisium = els.assay.value === "visium";
  els.visiumOptions.hidden = !isVisium;
  els.flexOptions.hidden = isVisium;

  if (!isVisium) {
    const flexV1 = els.flexVersion.value === "v1";
    if (!flexV1) {
      els.flexMode.value = "singleplex";
    }
    els.flexMode.disabled = !flexV1;
    const multiplex = flexV1 && els.flexMode.value === "multiplex";
    els.flexBarcode.disabled = !multiplex;
    els.flexNN.disabled = !flexV1;
  }
}

function refreshTable() {
  const filtered = filterRows(lastDesign.rows, els.search.value);
  const sorted = sortRows(filtered, sortState.col, sortState.dir);
  currentViewRows = sorted;
  renderTable(els.table, currentViewRows, TABLE_COLUMNS, sortState);
}

els.table.addEventListener("click", (e) => {
  const th = e.target.closest("th");
  if (!th) return;
  const col = th.dataset.col;
  if (!col) return;
  if (sortState.col === col) {
    sortState.dir = sortState.dir === "asc" ? "desc" : "asc";
  } else {
    sortState.col = col;
    sortState.dir = "asc";
  }
  refreshTable();
});

els.search.addEventListener("input", () => refreshTable());
els.assay.addEventListener("change", () => updateAssayControls());
els.flexVersion.addEventListener("change", () => updateAssayControls());
els.flexMode.addEventListener("change", () => updateAssayControls());

els.example.addEventListener("click", () => {
  els.fasta.value = `>target1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>target2
TTGACCGTATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA`;
});

els.file.addEventListener("change", async () => {
  const f = els.file.files?.[0];
  if (!f) return;
  const txt = await f.text();
  els.fasta.value = txt;
});

els.run.addEventListener("click", () => {
  try {
    const statusEl = els.status ?? document.getElementById("status");
    const statsEl = els.stats ?? document.getElementById("stats");
    if (statusEl) statusEl.textContent = "Running...";
    const params = getParams();
    const { rows, stats } = designProbes(els.fasta.value, params);
    lastDesign = { rows, stats };

    if (statsEl) {
      statsEl.textContent =
        `targets: ${stats.targets}\n` +
        `total_bases: ${stats.total_bases}\n` +
        `candidates (post-filters): ${stats.candidates}\n` +
        `selected: ${stats.selected}\n`;
    }

    els.dlCsv.disabled = rows.length === 0;
    els.dlOpools.disabled = rows.length === 0;

    sortState = { col: "probe_id", dir: "asc" };
    refreshTable();
    if (statusEl) {
      statusEl.textContent = rows.length ? "Done." : "No probes passed filters.";
    }
  } catch (err) {
    console.error(err);
    const statusEl = els.status ?? document.getElementById("status");
    if (statusEl) {
      statusEl.textContent = err?.message ? `Error: ${err.message}` : "Error (see console).";
    }
  }
});

els.dlCsv.addEventListener("click", () => {
  const csv = toCsv(lastDesign.rows, TABLE_COLUMNS);
  downloadText("probe_table.csv", csv);
});

els.dlOpools.addEventListener("click", () => {
  const poolName = els.poolName.value.trim() || "poolOne";
  const opRows = toOpoolsRows(lastDesign.rows, poolName);
  const cols = Object.keys(opRows[0] ?? { "Pool name": "", Sequence: "" });
  const csv = toCsv(opRows, cols);
  downloadText("opools.csv", csv);
});

// initial render
renderTable(els.table, [], TABLE_COLUMNS, sortState);
updateAssayControls();
