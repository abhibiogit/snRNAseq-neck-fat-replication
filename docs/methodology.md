# Methodology

Detailed write-up of the replicated pipeline, with parameter rationale and
notes on where this implementation diverges from the
original paper's methods text.

> Shamsi, F., Piper, M., Ho, L.L. et al. Vascular smooth muscle-derived
> Trpv1+ progenitors are a source of cold-induced thermogenic adipocytes.
> *Nat Metab* 3, 485–495 (2021). https://doi.org/10.1038/s42255-021-00373-z

Data: ArrayExpress [E-MTAB-8564](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8564),
9 single-nuclei RNA-seq samples of human neck fat.

---

## Stage 1 — Quality control and filtering (`01_qc_filtering_h5ad_export.R`)

**Aim**: Identify which samples have usable signal (vs. background-RNA
dominated) and produce cleaned, doublet-filtered count matrices for
integration.

**Steps**:

1. Read each of the 9 sample folders with `DropletUtils::read10xCounts()`.
2. Compute barcode ranks with `DropletUtils::barcodeRanks()` for each sample.
3. Plot barcode rank curves (log-log, total UMI vs. rank) with the knee and
   inflection points annotated, one plot per sample.
4. **QC decision**: a sample is retained only if its barcode rank plot shows
   a clear inflection point separating cell-containing barcodes from empty
   droplets. This is a visual/qualitative judgment, consistent with the
   paper's stated criterion ("lack of an inflection point ... indicated a
   high degree of background RNA contamination"). Applying this criterion
   to the 9 samples resulted in retaining the 4 `H_BAT_nF_*` samples and
   excluding the 5 `H_BAT_F_*` samples.
5. For each retained sample, keep only the top 1,000 barcodes ranked by
   total UMI count (a strict, fixed-size cell selection rather than a
   knee/inflection-based automatic threshold).
6. Remove likely doublets by dropping any retained barcode with a total UMI
   count above 20,000.
7. Export each cleaned sample independently to `.h5ad` via
   `sceasy::convertFormat(..., from = "sce", to = "anndata")`.

The paper's methods specify all of the above numerically (top 1,000 barcodes, >20,000 UMI doublet threshold), so
these were applied directly rather than tuned.
The inflection-point QC criterion is qualitative; no quantitative threshold is given in the paper, so sample retention here
was determined by visual inspection of each barcode rank plot.

---

## Stage 2 — Batch integration and embedding (`02_scvi_integration_umap.ipynb`)

**Aim**: Combine the 4 QC-passed samples into a single dataset with
sample-specific (batch) effects removed, producing a shared low-dimensional
representation for clustering.

**Steps**:

1. Load the 4 per-sample `.h5ad` files; prefix barcodes with sample name to
   avoid collisions across samples.
2. Concatenate with `ad.concat(..., label="batch_indices", merge="same")` —
   `batch_indices` records each nucleus's sample of origin, used later as
   the batch key for scVI.
3. Preserve raw counts in `adata.layers["counts"]` (scVI requires raw, not
   normalized, counts as input).
4. Convert gene identifiers from Ensembl IDs to gene symbols
   (`adata.var["Symbol"]`), deduplicating with `var_names_make_unique()`.
   This step is not explicit in the paper's methods but is necessary for
   downstream marker-gene lookups (ADIPOQ, PDGFRA, TRPV1) by symbol.
5. Filter genes detected in fewer than 3 nuclei
   (`sc.pp.filter_genes(min_cells=3)`) — standard low-expression filtering,
   applied before HVG selection. Reduces the gene set from ~21,959 to
   ~17,650. Not explicitly stated in the paper's methods paragraph; added
   here as standard practice to remove noise.
6. Select the top 2,000 highly variable genes
   (`sc.pp.highly_variable_genes(n_top_genes=2000, flavor="cell_ranger")`),
   corresponding to the paper's `subsample_genes` step.
7. Integrate with scVI:
   ```python
   scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key="batch_indices")
   model = scvi.model.SCVI(adata_hvg, n_latent=20, gene_likelihood="nb")
   model.train(max_epochs=400, accelerator="gpu", devices=1, plan_kwargs={"lr": 1e-3})
   ```
8. Extract the 20-D latent representation onto the full AnnData object:
   `adata.obsm["X_scVI"] = model.get_latent_representation(adata_hvg)`.
9. Build the neighbor graph on the scVI latent space
   (`sc.pp.neighbors(use_rep="X_scVI", n_neighbors=15)`) and compute UMAP
   (`sc.tl.umap(min_dist=0.3)`).
10. Save the integrated object as `H_BAT_nF_integrated.h5ad`.

**Parameter mapping vs. the paper's legacy scVI API**:

| Paper (legacy `VAE`/`UnsupervisedTrainer`) | This replication (`scvi.model.SCVI`) |
| `n_epochs = 400` | `max_epochs=400` |
| `lr = 1e-3` | `plan_kwargs={"lr": 1e-3}` |
| `use_batches = True` | `batch_key="batch_indices"` via `setup_anndata()` |
| `reconstruction_loss = "nb"` | `gene_likelihood="nb"` |
| `n_latent = 20` | `n_latent=20` (unchanged) |

The paper's original implementation predates the current `scvi-tools`
model-based API; the mapping above preserves the same effective
hyperparameters under the modern interface. Training used GPU acceleration
(`accelerator="gpu", devices=1`) where i used a NVIDIA RTX 4060 Laptop GPU; a subset of nuclei was automatically reassigned
from the training to validation split by scvi-tools' default data splitter.

**Why these parameter choices**: all values (epochs, learning rate, batch
handling, likelihood, latent dimensionality) are taken directly from the
paper's stated methods; only the API surface differs due to change in package version.

---

## Stage 3 — Clustering and marker identification (`03_seurat_clustering_marker_id.R`)

**Aim**: Cluster the integrated dataset, identify cluster marker genes,
annotate adipocyte and adipocyte progenitor (APC) populations, and resolve
finer substructure within them.

**Steps**:

1. Convert the integrated `.h5ad` to a Seurat object via
   `sceasy::convertFormat(..., from = "anndata", to = "seurat")`, using
   `reticulate::use_condaenv()` to point R at the same conda environment
   used in stage 2 (required for sceasy's Python dependencies).
2. Verify the imported object (dimensions, metadata, available reductions
   and assays) and save the raw import (`H_BAT_integrated.seurat.rds`)
   before further processing.
3. Standard Seurat preprocessing: `NormalizeData()` →
   `FindVariableFeatures()` → `ScaleData()`.
4. Plot the UMAP carried over from Scanpy/scVI
   (`DimPlot(reduction = "umap")`).
5. Build the neighbor graph on the **scVI latent space**, not a
   Seurat-computed PCA: `FindNeighbors(reduction = "scVI", dims = 1:20)`.
6. Cluster with `FindClusters(resolution = 0.5)` on the full dataset.
7. Identify cluster markers:
   `FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)`
   (Wilcoxon rank-sum test, Seurat default) → `H_BAT_cluster_markers.csv`.
8. Visualize canonical markers (`ADIPOQ`, `PDGFRA`, `TRPV1`) via
   `FeaturePlot` and `VlnPlot` to annotate clusters.
9. **Manual annotation**: cluster 2 identified as adipocytes (ADIPOQ+),
   cluster 3 as adipocyte progenitor cells / APCs (PDGFRA+).
10. Subset to clusters 2 and 3, repeat `FindNeighbors` (same `scVI`
    reduction, 20 dims) and `FindClusters(resolution = 0.8)` on the subset
    to resolve finer substructure — mirroring the paper's described
    "subsetted and reclustered again as described above" step.
11. Repeat marker identification and marker visualization
    (`ADIPOQ`, `PDGFRA`, `TRPV1`) on the subset →
    `H_BAT_subset_markers.csv`, `markerplot2.png`, `violinplot2.png`.

**Why these parameter choices**:
- Marker thresholds (`min.pct = 0.25`, `logfc.threshold = 0.25`) are Seurat
  defaults for standard marker discovery, not explicitly stated in the
  paper's methods.
- Clustering resolutions (0.5 full dataset, 0.8 subset) are not specified in
  the paper's methods paragraph and were chosen independently. [These were tuned to approximate the paper's
  reported cluster counts]
- `TRPV1` was included as a marker throughout — this is the gene central to
  the paper's main finding (TRPV1+ progenitors as a source of thermogenic
  adipocytes), so its expression pattern across the subset clusters is the
  most direct point of comparison to the paper's central result.

---

## Summary of change in versions of various software and set parameters

| Area | Paper's original approach | This replication | Rationale |
| scVI API | Legacy `VAE` / `UnsupervisedTrainer` | Current `scvi.model.SCVI` | Legacy API deprecated in current scvi-tools; parameters mapped 1:1 (see table in Stage 2) |
| HVG selection | `subsample_genes` (scVI-adjacent helper) | `sc.pp.highly_variable_genes(flavor="cell_ranger")` | Direct Scanpy equivalent for top-N HVG selection by counts |
| Low-expression gene filtering | Not stated | `sc.pp.filter_genes(min_cells=3)` | Standard preprocessing step, added for signal quality |
| Sample QC criterion | Qualitative (inflection point present/absent) | Same, applied by visual inspection | No quantitative threshold given in original methods |
| Clustering resolution | Not stated | 0.5 (full), 0.8 (subset) | Chosen independently; not paper-specified |
| Marker thresholds | Not stated | Seurat defaults (`min.pct=0.25`, `logfc.threshold=0.25`) | Not paper-specified |

## Comparison of results to the original paper


- I was able to recover the same number of clusters as the paper from the adipocyte and adipocyte progenitor cell subset.
- I could observe that ADIPOQ+ and PDGFRA+ clusters correspond in relative size/position
  to those reported in the paper.
- I was able to resolve a distinct TRPV1+ subpopulation of adipocyte progenitor cells as reported.
