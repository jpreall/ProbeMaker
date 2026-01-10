# ProbeMaker

Browser-based designer for 10x Genomics Visium and Chromium Flex spike-in probes. Paste or upload FASTA, set assay options, and generate order-ready LHS/RHS oligos with required constant arms. No BLAST/off-target screening is performed in this web version.

## Features
- 50 nt split probes (25/25) with GC and homopolymer constraints
- Ligation constraint: LHS probe 3' base = T (optional)
- Visium and Flex (v1, v2, v2 4-plex) constant arms appended
- CSV export of full probe table
- IDT oPools export with LHS first, then RHS

## Quick start
1. Open `index.html` in a browser (or use the GitHub Pages URL).
2. Paste FASTA or upload a `.fa/.fasta/.fna/.txt` file.
3. Choose assay options and constraints.
4. Click "Design probes".
5. Download the table CSV or oPools CSV.

## Inputs and constraints
- **FASTA**: one or more sequences (A/C/G/T/N, U allowed and treated as T).
- **GC min/max**: applied separately to LHS and RHS 25-mers.
- **Max homopolymer run**: applied across the concatenated LHS+RHS.
- **Min spacing**: minimum distance between probe start positions.
- **Auto-pick count**: number of probes selected per target. Blank = all candidates.
- **Ligation constraint**: requires LHS probe to end with T at the 3' junction.

## Assay options
- **Visium**: single option (v1/v2/HD) with shared design and constants.
- **Flex v1**: singleplex or multiplex; multiplex requires BC001..BC016 and NN.
- **Flex v2**: uses GEM-X constants; no barcode/NN.
- **Flex v2 4-plex**: uses pCS1 RHS constant.

## File output
- `probe_table.csv`: full table including sequences and order-ready oligos.
- `opools.csv`: two rows per probe (`LHS` then `RHS`) with a single pool name column.

## Notes
- This web tool is meant for rapid probe design and ordering.
- For BLAST/off-target filtering, use `ProbeDesign10X.py`.

## License
No license specified. Add one if you intend to distribute this broadly.
