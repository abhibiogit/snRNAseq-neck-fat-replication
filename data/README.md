# Data

Raw data is **not** committed to this repository (public data, large file
sizes).

## Source

- **Accession**: E-MTAB-8564
- **Repository**: [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8564)
- **Description**: Single-nuclei RNA-seq of human neck (brown/beige) fat, 9
  samples, used in:
  Shamsi, F., Piper, M., Ho, L.L. et al. Vascular smooth muscle-derived
  Trpv1+ progenitors are a source of cold-induced thermogenic adipocytes.
  *Nat Metab* 3, 485–495 (2021). https://doi.org/10.1038/s42255-021-00373-z

## How to download

1. Go to the ArrayExpress accession page above.
2. Download the CellRanger output folders for all 9 samples (or use the FTP
   links provided on the accession page).
3. Place downloaded sample folders under `data/raw/` (this path is
   gitignored), matching the folder names expected by
   `scripts/01_qc_filtering_h5ad_export.R`:

```bash
data/raw/
├── H_BAT_F_1/
├── H_BAT_F_5/
├── H_BAT_F_6/
├── H_BAT_F_7/
├── H_BAT_F_8/
├── H_BAT_nF_1/
├── H_BAT_nF_2/
├── H_BAT_nF_3/
└── H_BAT_nF_4/
```

Update the `base_dir` variable in `scripts/01_qc_filtering_h5ad_export.R` to
point at `data/raw/` (or wherever you place these).

## Samples used in this replication

Of the 9 samples, the following 4 passed the barcode-rank QC check
(see `scripts/01_qc_filtering_h5ad_export.R`) and were used downstream:

- H_BAT_nF_1
- H_BAT_nF_2
- H_BAT_nF_3
- H_BAT_nF_4

The remaining 5 (`H_BAT_F_1`, `H_BAT_F_5`, `H_BAT_F_6`, `H_BAT_F_7`,
`H_BAT_F_8`) were QC'd and plotted but excluded from downstream analysis due
to a lack of a clear inflection point in the barcode rank plot (high
background RNA contamination), consistent with the original paper's
exclusion criteria.
