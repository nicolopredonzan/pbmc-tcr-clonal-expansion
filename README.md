# T-cell clonal expansion in a healthy donor (10x 5' GEX + V(D)J, PBMC)

Small exploratory re-analysis of a public 10x Genomics dataset, written to practise paired
single-cell RNA + TCR analysis in R/Seurat. It is **not** a new finding: one donor, descriptive only.

## Question
In a healthy donor, which T-cell states carry expanded TCR clones (identical CDR3α + CDR3β in ≥3 cells)?
Expectation, fixed before looking at the clone counts: expansions sit in CD8 memory/effector cells, not in naive T cells.

## Data
- 10x Genomics, **"10k Human PBMCs, 5' v2.0, Chromium X (without intronic reads)"**, CC BY 4.0
  ([dataset page](https://www.10xgenomics.com/jp/datasets/10k-human-pbmcs-5-v2-0-chromium-x-without-intronic-reads-2-standard)).
- Healthy donor (PBMC from AllCells), 5' gene expression + V(D)J T from the same cells, Cell Ranger 6.1.2 (`cellranger multi`).
- Files used (not redistributed, download them into `data/`):
  `..._count_sample_feature_bc_matrix.h5` and `..._vdj_t_filtered_contig_annotations.csv`
  (exact file names in the config block of the script).

## Run
```r
source("pbmc_tcr_clonal_expansion.R")   # Seurat v5, hdf5r, dplyr, ggplot2 (+ ggrepel, patchwork, scales, installed with Seurat)
```
Optional (recommended): `BiocManager::install("scDblFinder")` for doublet flagging; the script skips this step if the package is missing.

Outputs: `results/figures/` (PNG, 300 dpi), `results/tables/` (CSV tables, labelled Seurat object, `sessionInfo.txt`). Seed 42.
Tested with R 4.5.2 and Seurat 6.0.0 (v5 API); exact versions in `results/tables/sessionInfo.txt`.

## Methods in brief
- **QC**: 500 < genes < 4000 per cell, mitochondrial reads < 15%. Thresholds were chosen from violin/scatter plots and
  by counting how many cells with a complete TCR each threshold removes (`results/tables/qc_threshold_grid.csv`, heatmap
  `results/figures/04_qc_threshold_loss_heatmap`). The upper gene cap removes no cell belonging to an expanded clone.
  9,331 → 7,350 cells.
- **Doublets**: flagged with scDblFinder but not removed (removing cells would change the clustering and the
  hard-coded labels); used to judge the expanded cells that fall in non-T clusters (`results/tables/doublets_by_cell_type.csv`).
- **Clustering**: 2,000 variable genes (TCR genes excluded so clusters are not driven by clonotype), 15 PCs (elbow),
  resolution 0.3. Resolution 0.6 is kept as a comparison (`results/tables/resolution_comparison.csv`): it splits clusters further
  without being needed for this question.
- **Labels**: manual, from canonical markers (DotPlot, a second DotPlot for minor populations, FindMarkers).
  Minor populations: plasmacytoid DC (SCT, LRRC26, LILRA4, GZMB, CD4, no CD3E) and platelets (GP9, PPBP, PF4). The cluster-to-label mapping is hard-coded and checked
  by `check_marker()`, because cluster numbers are not stable across versions/seeds. Naive T cells cannot be split into
  CD4/CD8 by RNA alone (no surface-protein data in the h5).
- **Clone**: identical CDR3α + CDR3β, only in cells with exactly one productive TRA and one TRB contig (3,481 cells after QC).
  Cells with two alpha chains or a missing chain are not assigned to a clone; cells with two TRB are probable doublets.
  **Expanded** = clone of ≥3 cells (size counted on all cells after QC). Clones are defined on CDR3 amino acids only;
  adding identical TRAV/TRBV genes is reported as a sensitivity check (`results/tables/clone_definition_comparison.csv`).
- **MAIT check**: per cluster, fraction of TCR 1+1 cells with the semi-invariant MAIT alpha chain (TRAV1-2 with TRAJ33/12/20)
  and of KLRB1/SLC4A10/ZBTB16-positive cells (`results/tables/mait_check_by_cluster.csv`). Cluster 9 ("GZMK+ T": CD8 low,
  IL7R+, GZMK+) is the candidate; it is excluded in the sensitivity definition of the CD8 memory/effector group.
- **Repertoire**: clone size distribution (`18_clone_size_distribution`) and Shannon diversity per T-cell type,
  also after downsampling to the smallest group (`results/tables/repertoire_diversity_by_cell_type.csv`).

## Results (descriptive)
- 29 expanded clones; the largest has 14 cells. 126 of 3,481 cells (3.6%) belong to an expanded clone.
- 110/534 cells (21%) in the CD8 GZMK+ memory/effector clusters belong to an expanded clone
  (104/463, 22%, if the second GZMK+ cluster is dropped), versus **0/1,153 in the naive T clusters**.
  No expanded clone has any cell in a naive cluster.
- CD4 T cluster (mixed naive/memory): 11/1,614 cells (0.7%).
- 16 cells of expanded clones sit outside these groups: 11 in the CD4 T cluster, 2 in CD14+ monocytes, 1 NK, 1 pDC,
  1 unassigned (`results/tables/expanded_cells_outside_test_groups.csv`, not removed). Three small clones (3 cells each)
  sit entirely in the CD4 T cluster and are plausible CD4 memory expansions. The single cells with monocyte/NK/pDC labels
  belong to clones whose other cells are CD8 memory/effector: most plausibly doublets or ambient TCR transcripts
  (see the `dbl_class` column).

- Sanity check independent of the cluster labels (`results/figures/16_dotplot_expanded_vs_not`): among TCR 1+1 cells, those in
  expanded clones show higher CCL5/NKG7/GZMK/GZMH/CD8A and lower CCR7/SELL/TCF7 than non-expanded cells.

Key figures:

![Cell types and expanded clones (all cells)](results/figures/14_umap_expanded_clones.png)

![Summary: expanded cells on the UMAP and expansion by cell type](results/figures/15_summary_expanded_clones.png)

Zoom on the T-cell clusters only (non-T clusters hidden; the few cells of expanded clones in those clusters are in `results/tables/expanded_cells_outside_test_groups.csv`):

![T-cell clusters and expanded clones](results/figures/17_umap_T_cells_expanded_clones.png)

## Sample identity check: donor sex
Checking sex-linked genes is a cheap, standard control that the cells are the sample the metadata describe
(in multi-donor studies it is also used to catch sample swaps or to split pooled donors).

- **Metadata**: the 10x Genomics dataset page describes the donor as **female**.
- **Data**: XIST (expressed from the inactive X, so present only in female cells) is detected in **0%** of cells
  after QC, while Y-linked genes are detected in most cells: RPS4Y1 88%, DDX3Y 34%, USP9Y 24%
  (`results/tables/sex_check.csv`). This pattern is that of **male** cells, not of a female sample with noise.
- **Interpretation**: I cannot tell from the data alone whether the metadata on the dataset page are wrong or the
  files come from a different donor than described. I report the discrepancy and do not correct it.
- **Impact on this analysis**: none for the question asked (expansion vs. cell state *within* one donor; no sex-based
  comparison is made). It would matter if this sample were combined with others or compared by sex.

The general lesson I took from it: public metadata should be verified against the data before being used.

## Limitations
- One donor, no replicate: no inferential statistics, no claim of generality.
- The result is the expected one (memory/effector expansion): the project shows method, not discovery.
- Marker-based labels; the CD4 T cluster mixes naive and memory cells.
- Cells without a complete TCR 1+1 are not analysed; TCR capture may differ between cell states.
- Donor sex does not match the metadata (see "Sample identity check" below).
- Clones are defined on CDR3 amino acid sequences; identical CDR3s from different V genes or nucleotide sequences
  (convergent recombination, e.g. in MAIT cells) are counted as one clone.

## Use of AI
The code was developed with the help of an AI assistant (Claude). I ran the analysis, chose the parameters
and can explain each step.
