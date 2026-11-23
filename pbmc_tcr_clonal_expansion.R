# =============================================================================
# Clonal expansion vs transcriptional state in healthy-donor PBMC
# Paired 5' scRNA-seq + scTCR-seq (10x Genomics), R / Seurat (v5 API, run with Seurat 6.0.0)
# =============================================================================
# Question (fixed before looking at clone counts per group):
#   Are expanded T-cell clones (>= 3 cells sharing the same TCR) found in
#   CD8 memory/effector T cells and not in naive T cells?
#
# Data: 10x Genomics public dataset (CC BY 4.0)
#   "10k Human PBMCs, 5' v2.0, Chromium X (without intronic reads)"
#   https://www.10xgenomics.com/jp/datasets/10k-human-pbmcs-5-v2-0-chromium-x-without-intronic-reads-2-standard
#   Healthy donor, Cell Ranger 6.1.2 (cellranger multi). Not redistributed here.
#   Put the two files below in ./data/
#     - ..._count_sample_feature_bc_matrix.h5
#     - ..._vdj_t_filtered_contig_annotations.csv
#
# Definitions
#   - Clone: identical CDR3-alpha + CDR3-beta pair (amino acids), only in cells
#     with exactly one productive TRA and one productive TRB contig (TCR "1+1").
#     Cells with two alpha chains or a missing chain are NOT assigned to a clone.
#     Sensitivity check: same definition plus identical TRAV and TRBV genes.
#   - Expanded clone: >= EXPANDED_MIN_SIZE cells (counted after QC, all cells).
#
# Limitations (state them in the README)
#   - One donor; descriptive analysis only, no p-values are interpreted
#     (cells are not independent samples).
#   - Cell types are annotated from marker genes only (RNA, no surface protein).
#   - The "Naive T" cluster mixes CD4 and CD8 naive cells (not resolved by RNA).
#   - Cluster-to-label mapping is hard-coded for this exact run (seed, Seurat
#     version, parameters); sanity checks below stop the script if it breaks.
#   - Doublets are flagged (scDblFinder) but NOT removed: removing cells would
#     change the clustering and therefore the hard-coded cluster labels.
#   - Donor sex: the 10x page reports a female donor, but XIST is not detected
#     and Y-linked genes are detected in most cells (results/tables/sex_check.csv):
#     the expression data are male. Metadata error or different donor cannot be
#     distinguished from the data; reported, not corrected. No effect on this
#     within-donor question.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(hdf5r)
})
set.seed(42)

# ---- Config -----------------------------------------------------------------
H5_PATH     <- "data/10k_PBMC_5pv2_nextgem_Chromium_X_10k_PBMC_5pv2_nextgem_Chromium_X_count_sample_feature_bc_matrix.h5"
CONTIG_PATH <- "data/10k_PBMC_5pv2_nextgem_Chromium_X_10k_PBMC_5pv2_nextgem_Chromium_X_vdj_t_filtered_contig_annotations.csv"
FIG_DIR <- file.path("results", "figures")
RES_DIR <- file.path("results", "tables")

QC_MIN_FEATURES <- 500    # genes per cell
QC_MAX_FEATURES <- 4000
QC_MAX_MT       <- 15     # % mitochondrial reads
N_HVG           <- 2000
N_PCS           <- 15     # from the elbow plot
RESOLUTION      <- 0.3    # final clustering resolution
RESOLUTION_ALT  <- 0.6    # only for comparison
EXPANDED_MIN_SIZE <- 3

# Test groups are defined from cluster IDs at RESOLUTION (chosen from markers)
NAIVE_CLUSTERS           <- c("2", "5")
EFFECTOR_CLUSTERS        <- c("3", "9")   # main definition
EFFECTOR_CLUSTERS_STRICT <- c("3")        # sensitivity: cluster 9 may be MAIT cells (see MAIT check)

# MAIT cells: semi-invariant TCR alpha chain TRAV1-2 joined to TRAJ33 (canonical),
# TRAJ12 or TRAJ20
MAIT_TRAV <- "TRAV1-2"
MAIT_TRAJ <- c("TRAJ33", "TRAJ12", "TRAJ20")

RUN_MARKERS  <- TRUE   # FindMarkers for uncertain clusters (a few minutes)
RUN_DOUBLETS <- TRUE   # scDblFinder (Bioconductor); skipped if not installed

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, width = 7, height = 5) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
}

# ---- 1. Load RNA ------------------------------------------------------------
gex <- Read10X_h5(H5_PATH)
if (is.list(gex)) gex <- gex[["Gene Expression"]]   # if several modalities
obj <- CreateSeuratObject(gex, min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
message("Cells after loading: ", ncol(obj))

# ---- 2. Load TCR contigs ----------------------------------------------------
ctg <- read.csv(CONTIG_PATH, stringsAsFactors = FALSE)
stopifnot(all(c("barcode", "chain", "productive", "cdr3", "v_gene", "j_gene") %in% colnames(ctg)))
print(table(chain = ctg$chain, productive = ctg$productive))

# note: 'productive' is lowercase "true" in this file
ctg <- ctg %>%
  filter(chain %in% c("TRA", "TRB"),
         tolower(as.character(productive)) == "true",
         cdr3 != "None", !is.na(cdr3))
message("Productive TRA/TRB contigs with CDR3: ", nrow(ctg))

# one row per cell; n_TRA / n_TRB count CONTIGS (not distinct CDR3 sequences)
per_cell <- ctg %>%
  group_by(barcode) %>%
  summarise(
    n_TRA = sum(chain == "TRA"),
    n_TRB = sum(chain == "TRB"),
    TRA   = paste(sort(unique(cdr3[chain == "TRA"])), collapse = ";"),
    TRB   = paste(sort(unique(cdr3[chain == "TRB"])), collapse = ";"),
    TRAV  = paste(sort(unique(v_gene[chain == "TRA"])), collapse = ";"),
    TRAJ  = paste(sort(unique(j_gene[chain == "TRA"])), collapse = ";"),
    TRBV  = paste(sort(unique(v_gene[chain == "TRB"])), collapse = ";"),
    .groups = "drop"
  )
message("Cells with a TCR: ", nrow(per_cell))
print(table(n_TRA = per_cell$n_TRA, n_TRB = per_cell$n_TRB))

# barcodes must match between RNA and TCR
n_match <- sum(per_cell$barcode %in% colnames(obj))
message("TCR barcodes found in RNA object: ", n_match, " / ", nrow(per_cell))
stopifnot(n_match == nrow(per_cell))

# ---- 3. Clone definition (before QC, used for the QC justification) --------
per_cell <- per_cell %>%
  mutate(paired_1_1  = n_TRA == 1 & n_TRB == 1,
         clone_id    = ifelse(paired_1_1, paste(TRA, TRB, sep = "|"), NA),
         clone_id_vj = ifelse(paired_1_1, paste(TRAV, TRA, TRBV, TRB, sep = "|"), NA))
md_tcr <- tibble::column_to_rownames(
  per_cell[, c("barcode", "n_TRA", "n_TRB", "paired_1_1", "clone_id", "clone_id_vj",
               "TRAV", "TRAJ", "TRBV")], "barcode")
obj <- AddMetaData(obj, md_tcr)
obj$has_TCR <- !is.na(obj$n_TRA)
obj$mait_tcr <- obj$paired_1_1 %in% TRUE & obj$TRAV %in% MAIT_TRAV & obj$TRAJ %in% MAIT_TRAJ
message("Cells with TCR 1+1: ", sum(obj$paired_1_1, na.rm = TRUE))

# ---- 4. QC plots (before filtering) ----------------------------------------
# (the Seurat warning about the "data" layer is harmless: not normalized yet)
save_plot(VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size = 0),
          "01_qc_violin_before", 10, 4.5)
save_plot(FeatureScatter(obj, "nFeature_RNA", "percent.mt"), "02a_qc_scatter_features_mt")
save_plot(FeatureScatter(obj, "nCount_RNA", "nFeature_RNA"), "02b_qc_scatter_counts_features")
save_plot(ggplot(obj@meta.data, aes(nFeature_RNA)) +
            geom_histogram(bins = 100) + theme_classic(),
          "03a_hist_nFeature")
save_plot(ggplot(obj@meta.data, aes(percent.mt)) +
            geom_histogram(bins = 100) + theme_classic(),
          "03b_hist_percent_mt")

# ---- 5. Justify QC thresholds: how many TCR 1+1 cells would be lost? -------
md <- obj@meta.data
grid <- expand.grid(max_mt   = c(8, 10, 12, 15, 20, 25, 30),
                    min_feat = c(200, 250, 500, 800, 1000))
removed <- function(mt, nf) md$percent.mt >= mt | md$nFeature_RNA < nf
grid$all_cells <- mapply(function(mt, nf) sum(removed(mt, nf)),
                         grid$max_mt, grid$min_feat)
grid$tcr_cells <- mapply(function(mt, nf) sum(removed(mt, nf) & md$paired_1_1 %in% TRUE),
                         grid$max_mt, grid$min_feat)
write.csv(grid, file.path(RES_DIR, "qc_threshold_grid.csv"), row.names = FALSE)

long <- rbind(
  data.frame(grid[, c("max_mt", "min_feat")], metric = "All cells removed",      n = grid$all_cells),
  data.frame(grid[, c("max_mt", "min_feat")], metric = "TCR 1+1 cells removed",  n = grid$tcr_cells))
long <- long %>% group_by(metric) %>% mutate(rel = n / max(n)) %>% ungroup()
long$max_mt   <- factor(long$max_mt)
long$min_feat <- factor(long$min_feat)
chosen <- long %>% filter(as.character(max_mt) == as.character(QC_MAX_MT),
                          as.character(min_feat) == as.character(QC_MIN_FEATURES))
p_grid <- ggplot(long, aes(min_feat, max_mt, fill = rel)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  geom_tile(data = chosen, fill = NA, color = "red", linewidth = 1) +
  facet_wrap(~ metric) +
  scale_fill_gradient(low = "white", high = "steelblue", guide = "none") +
  labs(x = "Minimum genes per cell (cells below are removed)",
       y = "Maximum % mitochondrial reads (cells at or above are removed)",
       title = "Cells removed by QC thresholds (red = chosen)",
       subtitle = "Upper gene cap not included; fill is relative within each panel") +
  theme_minimal()
save_plot(p_grid, "04_qc_threshold_loss_heatmap", 10, 4.5)

# upper gene cap: do the removed TCR 1+1 cells belong to expanded clones?
cs_pre <- table(md$clone_id)
big <- md$nFeature_RNA >= QC_MAX_FEATURES & md$paired_1_1 %in% TRUE
big_sizes <- as.integer(cs_pre[as.character(md$clone_id[big])])
message("Cells with >= ", QC_MAX_FEATURES, " genes: ", sum(md$nFeature_RNA >= QC_MAX_FEATURES),
        " | TCR 1+1: ", sum(big),
        " | in clones >=2: ", sum(big_sizes >= 2),
        " | in clones >=", EXPANDED_MIN_SIZE, ": ", sum(big_sizes >= EXPANDED_MIN_SIZE))

# ---- 6. Apply QC ------------------------------------------------------------
n_before <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA > QC_MIN_FEATURES &
                nFeature_RNA < QC_MAX_FEATURES & percent.mt < QC_MAX_MT)
message("Cells after QC: ", ncol(obj), " (removed ", n_before - ncol(obj), ")")
message("TCR 1+1 cells kept: ", sum(obj$paired_1_1, na.rm = TRUE))
save_plot(VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size = 0),
          "05_qc_violin_after", 10, 4.5)

# clone sizes recomputed AFTER QC, on all cells (not on a subset of clusters)
cs <- table(obj$clone_id)
obj$clone_size <- as.integer(cs[as.character(obj$clone_id)])
obj$expanded   <- !is.na(obj$clone_size) & obj$clone_size >= EXPANDED_MIN_SIZE
print(table(clone_size = as.integer(cs)))
message("Expanded clones (>= ", EXPANDED_MIN_SIZE, " cells): ",
        sum(as.integer(cs) >= EXPANDED_MIN_SIZE))

# sensitivity: does adding the V genes to the clone definition change the result?
cs_vj <- table(obj$clone_id_vj)
clone_def_cmp <- data.frame(
  definition        = c("CDR3a + CDR3b (main)", "TRAV + CDR3a + TRBV + CDR3b"),
  clones            = c(length(cs), length(cs_vj)),
  expanded_clones   = c(sum(cs >= EXPANDED_MIN_SIZE), sum(cs_vj >= EXPANDED_MIN_SIZE)),
  cells_in_expanded = c(sum(cs[cs >= EXPANDED_MIN_SIZE]), sum(cs_vj[cs_vj >= EXPANDED_MIN_SIZE])))
print(clone_def_cmp)
write.csv(clone_def_cmp, file.path(RES_DIR, "clone_definition_comparison.csv"), row.names = FALSE)

# clone size distribution
size_df <- as.data.frame(table(clone_size = as.integer(cs)))
save_plot(ggplot(size_df, aes(clone_size, Freq)) +
            geom_col(fill = "steelblue") +
            geom_text(aes(label = Freq), vjust = -0.3, size = 3) +
            scale_y_log10() + theme_classic() +
            labs(x = "Clone size (cells)", y = "Number of clones (log scale)",
                 title = "Clone size distribution (TCR 1+1 cells after QC)"),
          "18_clone_size_distribution", 6, 4)

# ---- 6b. Doublet detection (flag only, cells are NOT removed) ---------------
obj$dbl_class <- NA_character_
obj$dbl_score <- NA_real_
if (RUN_DOUBLETS && requireNamespace("scDblFinder", quietly = TRUE)) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = LayerData(obj, assay = "RNA", layer = "counts")))
  set.seed(42)
  sce <- scDblFinder::scDblFinder(sce)
  obj$dbl_class <- as.character(sce$scDblFinder.class)
  obj$dbl_score <- sce$scDblFinder.score
  print(table(doublet = obj$dbl_class))
} else {
  message("scDblFinder not run (install with BiocManager::install('scDblFinder')).")
}

# ---- 7. Normalisation, HVGs, PCA -------------------------------------------
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = N_HVG, verbose = FALSE)
# remove TCR V/D/J/C genes from HVGs, otherwise cells cluster by clonotype
vf <- VariableFeatures(obj)
VariableFeatures(obj) <- vf[!grepl("^TR[ABGD][VDJC]", vf)]
message("HVGs used for PCA: ", length(VariableFeatures(obj)))
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
save_plot(ElbowPlot(obj, ndims = 30), "06_elbow")

# ---- 8. Clustering ----------------------------------------------------------
obj <- FindNeighbors(obj, dims = 1:N_PCS, verbose = FALSE)
obj <- FindClusters(obj, resolution = c(RESOLUTION, RESOLUTION_ALT), verbose = FALSE)
res_col     <- paste0(DefaultAssay(obj), "_snn_res.", RESOLUTION)
res_col_alt <- paste0(DefaultAssay(obj), "_snn_res.", RESOLUTION_ALT)
xtab <- table(res_low = obj@meta.data[[res_col]], res_high = obj@meta.data[[res_col_alt]])
print(xtab)
write.csv(as.data.frame.matrix(xtab), file.path(RES_DIR, "resolution_comparison.csv"))

obj$cluster <- as.character(obj@meta.data[[res_col]])
Idents(obj) <- obj$cluster
obj <- RunUMAP(obj, dims = 1:N_PCS, verbose = FALSE)
save_plot(DimPlot(obj, group.by = "cluster", label = TRUE), "07_umap_clusters")

# ---- 9. Marker genes for annotation ----------------------------------------
marker_genes <- c("CD3E", "CD4", "CD8A", "CCR7", "SELL", "IL7R", "GZMK", "GZMH",
                  "GZMB", "NKG7", "FGFBP2", "MS4A1", "CD14", "LYZ")
save_plot(DotPlot(obj, features = marker_genes, group.by = "cluster") + RotatedAxis(),
          "08_dotplot_markers_by_cluster", 9, 6)
save_plot(FeaturePlot(obj, c("CD3E", "CD4", "CD8A", "CCR7", "GZMH", "NKG7")),
          "09_featureplot_markers", 10, 7)

# second panel for the minor / uncertain populations:
#   pDC (LILRA4, CLEC4C, IL3RA), platelets (PPBP, PF4, GP9), cDC (FCER1A, CD1C),
#   MAIT (KLRB1, SLC4A10, ZBTB16), gamma-delta T (TRDC, TRGC1/2),
#   CD8 TEMRA vs GZMK+ memory (GNLY, CX3CR1, FGFBP2), proliferating (MKI67)
extra_markers <- intersect(c("LILRA4", "CLEC4C", "IL3RA", "PPBP", "PF4", "GP9",
                             "FCER1A", "CD1C", "KLRB1", "SLC4A10", "ZBTB16",
                             "TRDC", "TRGC1", "TRGC2", "GNLY", "CX3CR1", "FGFBP2", "MKI67"),
                           rownames(obj))
save_plot(DotPlot(obj, features = extra_markers, group.by = "cluster") + RotatedAxis(),
          "08b_dotplot_extra_markers_by_cluster", 11, 6)

if (RUN_MARKERS) {
  # p-values are NOT interpreted (cells from one donor are not independent);
  # clusters are judged on log2FC and the fraction of expressing cells.
  for (cl_id in c("2", "4", "8", "9", "10", "11")) {
    m <- FindMarkers(obj, ident.1 = cl_id, only.pos = TRUE, min.pct = 0.25)
    m <- m[order(-m$avg_log2FC), ]
    write.csv(head(m, 20), file.path(RES_DIR, paste0("markers_cluster_", cl_id, ".csv")))
    message("Top markers, cluster ", cl_id, ":"); print(head(m[, c("avg_log2FC", "pct.1", "pct.2")], 10))
  }
}

# ---- 10. Sanity checks: donor sex and TCR coverage per cluster -------------
sex_genes <- intersect(c("XIST", "RPS4Y1", "DDX3Y", "USP9Y"), rownames(gex))
sex_frac  <- Matrix::rowMeans(gex[sex_genes, colnames(obj), drop = FALSE] > 0)
print(round(sex_frac, 3))
write.csv(data.frame(gene = names(sex_frac), fraction_cells_detected = as.numeric(sex_frac)),
          file.path(RES_DIR, "sex_check.csv"), row.names = FALSE)
message("XIST detected in ", sum(gex["XIST", ] > 0), " / ", ncol(gex), " cells (raw matrix)")

print(round(prop.table(table(cluster = obj$cluster, has_TCR = obj$has_TCR), 1), 2))

# ---- 11. Cell-type labels ---------------------------------------------------
# Based on marker genes (DotPlot, extra-marker DotPlot, FindMarkers).
#   cluster 10: SCT, LRRC26, GZMB, CD4, no CD3E  -> plasmacytoid DC
#   cluster 11: GP9, PVALB, TMEM40, CABP5        -> platelets
#   cluster 9 : CD8 low, IL7R, GZMK, NKG7        -> possibly MAIT (see MAIT check)
#   cluster 4 : CD3E+, IL7R+, no clear markers   -> left unassigned
lab <- c("0"  = "CD4 T (memory/naive)",
         "1"  = "CD14+ monocytes",
         "2"  = "Naive T (CD4/CD8 unresolved)",
         "3"  = "CD8 GZMK+ memory/effector",
         "4"  = "Unassigned",
         "5"  = "CD8 naive",
         "6"  = "B cells",
         "7"  = "NK cells",
         "8"  = "CD16+ monocytes",
         "9"  = "GZMK+ T",
         "10" = "pDC",
         "11" = "Platelets")
stopifnot(all(obj$cluster %in% names(lab)))
obj$cell_type <- factor(unname(lab[obj$cluster]), levels = unname(lab))

# the cluster numbers are not stable across versions/seeds: verify the labels
check_marker <- function(label, gene) {
  d <- FetchData(obj, c(gene, "cell_type"))
  m <- tapply(d[[gene]], d$cell_type, mean)
  if (names(which.max(m)) != label)
    stop("Label check failed: '", label, "' does not have the highest ", gene,
         ". Cluster numbering changed; re-annotate.")
}
check_marker("B cells", "MS4A1")
check_marker("CD14+ monocytes", "CD14")
check_marker("NK cells", "FGFBP2")
check_marker("pDC", "LILRA4")
check_marker("Platelets", "PPBP")

save_plot(DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE, label.size = 3),
          "10_umap_cell_types", 10, 6)
save_plot(DotPlot(obj, features = marker_genes, group.by = "cell_type") + RotatedAxis(),
          "11_dotplot_markers_by_label", 10, 6)

# test groups (fixed from markers, before looking at clone counts per group)
obj$test_group <- case_when(
  obj$cluster %in% NAIVE_CLUSTERS    ~ "Naive T",
  obj$cluster %in% EFFECTOR_CLUSTERS ~ "CD8 memory/effector T",
  TRUE ~ NA_character_)
obj$test_group_strict <- case_when(
  obj$cluster %in% NAIVE_CLUSTERS           ~ "Naive T",
  obj$cluster %in% EFFECTOR_CLUSTERS_STRICT ~ "CD8 memory/effector T",
  TRUE ~ NA_character_)
save_plot(DimPlot(obj, group.by = "test_group"), "12_umap_test_groups")

t_types <- c("CD4 T (memory/naive)", "Naive T (CD4/CD8 unresolved)", "CD8 naive",
             "CD8 GZMK+ memory/effector", "GZMK+ T")

# ---- 11b. MAIT check (TCR + RNA, per cluster) -------------------------------
# MAIT cells: TRAV1-2 with TRAJ33/12/20, KLRB1+ SLC4A10+ ZBTB16+.
mait_genes <- intersect(c("KLRB1", "SLC4A10", "ZBTB16"), rownames(obj))
mait_df <- cbind(obj@meta.data[, c("cluster", "cell_type", "paired_1_1", "TRAV", "mait_tcr")],
                 FetchData(obj, mait_genes))
mait_check <- mait_df %>%
  group_by(cluster, cell_type) %>%
  summarise(n_cells       = n(),
            n_TCR11       = sum(paired_1_1 %in% TRUE),
            frac_TRAV1_2  = mean(TRAV[paired_1_1 %in% TRUE] %in% MAIT_TRAV),
            frac_MAIT_TCR = mean(mait_tcr[paired_1_1 %in% TRUE]),
            across(all_of(mait_genes), ~ mean(.x > 0), .names = "frac_{.col}_pos"),
            .groups = "drop") %>%
  arrange(desc(frac_MAIT_TCR))
print(mait_check, width = Inf)
write.csv(mait_check, file.path(RES_DIR, "mait_check_by_cluster.csv"), row.names = FALSE)

# ---- 12. Clonal expansion by group -----------------------------------------
summarise_groups <- function(group_col) {
  obj@meta.data %>%
    filter(!is.na(.data[[group_col]]), paired_1_1 %in% TRUE) %>%
    group_by(group = .data[[group_col]]) %>%
    summarise(cells_TCR11       = n(),
              clones            = n_distinct(clone_id),
              expanded_clones   = n_distinct(clone_id[expanded]),
              cells_in_expanded = sum(expanded),
              .groups = "drop")
}
res_main   <- summarise_groups("test_group")
res_strict <- summarise_groups("test_group_strict")
cat("\n== Main definition (clusters", EFFECTOR_CLUSTERS, ") ==\n");        print(res_main)
cat("\n== Sensitivity (clusters", EFFECTOR_CLUSTERS_STRICT, "only) ==\n");  print(res_strict)
write.csv(res_main,   file.path(RES_DIR, "clones_by_group_main.csv"),        row.names = FALSE)
write.csv(res_strict, file.path(RES_DIR, "clones_by_group_sensitivity.csv"), row.names = FALSE)

# expanded cells in ALL cell types (not only the test groups)
by_type <- obj@meta.data %>%
  filter(paired_1_1 %in% TRUE) %>%
  group_by(cell_type) %>%
  summarise(cells_TCR11 = n(), cells_in_expanded = sum(expanded),
            frac_expanded = mean(expanded), .groups = "drop")
print(by_type)
write.csv(by_type, file.path(RES_DIR, "expanded_by_cell_type.csv"), row.names = FALSE)

# repertoire diversity per T-cell type. Shannon entropy depends on the number of
# cells, so it is also computed after random downsampling to the smallest group.
shannon <- function(ids) { p <- as.numeric(table(ids)) / length(ids); -sum(p * log(p)) }
div_df <- obj@meta.data %>%
  filter(paired_1_1 %in% TRUE, cell_type %in% t_types) %>%
  mutate(cell_type = droplevels(cell_type))
n_min <- min(table(div_df$cell_type))
set.seed(42)
rep_div <- div_df %>%
  group_by(cell_type) %>%
  summarise(cells_TCR11      = n(),
            clones           = n_distinct(clone_id),
            shannon          = shannon(clone_id),
            shannon_rarefied = mean(replicate(200, shannon(sample(clone_id, n_min)))),
            .groups = "drop") %>%
  mutate(n_rarefied = n_min)
print(rep_div)
write.csv(rep_div, file.path(RES_DIR, "repertoire_diversity_by_cell_type.csv"), row.names = FALSE)

# clone-level table (with V/J genes, MAIT-like TCR and doublet flags)
clones <- obj@meta.data %>%
  filter(paired_1_1 %in% TRUE) %>%
  group_by(clone_id) %>%
  summarise(size = n(),
            n_cell_types = n_distinct(cell_type),
            dominant_cell_type = names(which.max(table(cell_type))),
            TRAV = first(TRAV), TRAJ = first(TRAJ), TRBV = first(TRBV),
            mait_tcr = first(mait_tcr),
            n_doublets = sum(dbl_class %in% "doublet"),
            .groups = "drop") %>%
  arrange(desc(size))
write.csv(clones, file.path(RES_DIR, "clone_table.csv"), row.names = FALSE)
cat("\nClones spanning >1 cell type (expanded only): ",
    sum(clones$n_cell_types > 1 & clones$size >= EXPANDED_MIN_SIZE), "\n")
print(head(clones %>% filter(size >= EXPANDED_MIN_SIZE), 10), width = Inf)

# doublet fraction per cell type
if (any(!is.na(obj$dbl_class))) {
  dbl_by_type <- obj@meta.data %>%
    group_by(cell_type) %>%
    summarise(cells = n(), doublets = sum(dbl_class %in% "doublet"),
              frac_doublets = mean(dbl_class %in% "doublet"),
              expanded_cells = sum(expanded),
              expanded_doublets = sum(expanded & dbl_class %in% "doublet"),
              .groups = "drop")
  print(dbl_by_type)
  write.csv(dbl_by_type, file.path(RES_DIR, "doublets_by_cell_type.csv"), row.names = FALSE)
}

# ---- 13. Figures: expanded clones -------------------------------------------
# (a) fraction of cells in expanded clones, by cell type
frac_df <- by_type %>% filter(cells_TCR11 >= 20)
save_plot(ggplot(frac_df, aes(reorder(cell_type, frac_expanded), frac_expanded)) +
            geom_col(fill = "steelblue") + coord_flip() + theme_classic() +
            labs(x = NULL, y = paste0("Fraction of TCR 1+1 cells in clones >= ", EXPANDED_MIN_SIZE),
                 title = "Clonal expansion by cell type (cell types with >= 20 TCR 1+1 cells)"),
          "13_expanded_fraction_by_cell_type", 9, 4.5)

# (b) UMAP layers shared by the UMAP figures
um <- as.data.frame(Embeddings(obj, "umap"))
colnames(um) <- c("UMAP_1", "UMAP_2")
um$clone_id  <- obj$clone_id
um$expanded  <- obj$expanded
um$cell_type <- obj$cell_type

# cell-type names written on the map (used by figure 15; only cell types with >= 150 cells)
lab_pos <- um %>%
  group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), n = n(), .groups = "drop") %>%
  filter(n >= 150) %>%
  mutate(label = vapply(cell_type, function(x) paste(strwrap(x, 16), collapse = "\n"), ""))
add_labels <- function() {
  ggrepel::geom_text_repel(data = lab_pos, aes(UMAP_1, UMAP_2, label = label),
                           size = 5, fontface = "bold", color = "black",
                           bg.color = "white", bg.r = 0.3, lineheight = 0.9,
                           min.segment.length = Inf, box.padding = 0.8,
                           max.overlaps = Inf, seed = 42)
}

# Cell-type colours (muted; unassigned in grey) for the legend panels
ct_levels <- sort(unique(as.character(um$cell_type)))
ct_pal <- setNames(scales::hue_pal(c = 60, l = 70)(length(ct_levels)), ct_levels)
ct_pal["Unassigned"] <- "grey80"
um$cell_type <- factor(um$cell_type, levels = ct_levels)

exp_clones <- clones %>% filter(size >= EXPANDED_MIN_SIZE)
top_n   <- min(6, nrow(exp_clones))
top_ids <- exp_clones$clone_id[seq_len(top_n)]
other_name <- "Other expanded clones"
# Colour-blind-safe (Okabe-Ito) colours for the largest clones; all other
# expanded clones share one dark grey, so the legend can be read completely.
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2")
pal <- c(setNames(okabe[seq_len(top_n)], top_ids), setNames("grey35", other_name))
clone_lab <- c(paste0(top_ids, " (n=", exp_clones$size[seq_len(top_n)], ")"),
               paste0(other_name, " (", nrow(exp_clones) - top_n, " clones)"))

um_e <- um[um$expanded, ]
um_e$grp <- factor(ifelse(um_e$clone_id %in% top_ids, um_e$clone_id, other_name),
                   levels = c(top_ids, other_name))
um_e <- um_e[order(um_e$grp, decreasing = TRUE), ]   # largest clones drawn on top

# helper: panel A = cell types (colour + legend, no text on the map)
plot_types <- function(d, title) {
  ggplot(d, aes(UMAP_1, UMAP_2, color = cell_type)) +
    geom_point(size = 0.4) +
    scale_color_manual(values = ct_pal, name = "Cell type", drop = TRUE) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme_classic() + labs(x = "umap_1", y = "umap_2", title = title)
}
# helper: panel B = clones (six largest coloured, others dark grey)
plot_clones <- function(d, de, title, caption = NULL) {
  ggplot() +
    geom_point(data = d[!d$expanded, ], aes(UMAP_1, UMAP_2), color = "grey88", size = 0.3) +
    geom_point(data = de, aes(UMAP_1, UMAP_2, fill = grp),
               shape = 21, color = "black", stroke = 0.25, size = 2) +
    scale_fill_manual(values = pal, labels = clone_lab, name = "Expanded clones\n(CDR3a|CDR3b)",
                      drop = FALSE) +
    theme_classic() + labs(x = "umap_1", y = "umap_2", title = title, caption = caption)
}

# Figure 14: all cells. A = cell types, B = expanded clones (same coordinates)
p14a <- plot_types(um, "Cell types")
p14b <- plot_clones(um, um_e, paste0(nrow(exp_clones), " expanded clones (>= ",
                                     EXPANDED_MIN_SIZE, " cells); six largest coloured"))
save_plot(patchwork::wrap_plots(p14a, p14b, widths = c(1, 1.25)) +
            patchwork::plot_annotation(tag_levels = "A"),
          "14_umap_expanded_clones", 18, 6)

# Figure 17: zoom on T-cell clusters only (the question of the project)
um_t   <- um[um$cell_type %in% t_types, ]
um_t$cell_type <- droplevels(um_t$cell_type)
um_te  <- um_e[um_e$cell_type %in% t_types, ]
n_out  <- sum(um$expanded & !um$cell_type %in% t_types)
p17a <- plot_types(um_t, "T-cell clusters")
p17b <- plot_clones(um_t, um_te, "Expanded clones within T-cell clusters",
                    caption = paste0("Non-T clusters not shown (", n_out,
                                     " cells of expanded clones, see results/tables)"))
save_plot(patchwork::wrap_plots(p17a, p17b, widths = c(1, 1.25)) +
            patchwork::plot_annotation(tag_levels = "A"),
          "17_umap_T_cells_expanded_clones", 18, 6)

# (c) Main summary figure: annotated map + cells in expanded clones + quantification
um_flag <- um[order(um$expanded), ]                  # expanded cells drawn on top
p_flag <- ggplot(um_flag, aes(UMAP_1, UMAP_2, color = expanded, size = expanded)) +
  geom_point() +
  scale_color_manual(values = c(`FALSE` = "grey88", `TRUE` = "firebrick"),
                     labels = c("Other cells", "Cell in an expanded clone"), name = NULL) +
  scale_size_manual(values = c(`FALSE` = 0.3, `TRUE` = 1.4), guide = "none") +
  add_labels() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() + theme(legend.position = "bottom") +
  labs(x = "umap_1", y = "umap_2", title = "Cells in expanded clones")
p_frac <- ggplot(frac_df, aes(reorder(cell_type, frac_expanded), frac_expanded)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_classic() +
  labs(x = NULL, y = paste0("Fraction of TCR 1+1 cells in clones >= ", EXPANDED_MIN_SIZE),
       title = "Expansion by cell type")
save_plot(patchwork::wrap_plots(p_flag, p_frac, widths = c(1.3, 1)) +
            patchwork::plot_annotation(tag_levels = "A"),
          "15_summary_expanded_clones", 13, 5.5)

# (d) Sanity check independent of the cluster labels: do cells of expanded clones
# look like effector/memory cells (vs. naive) among all TCR 1+1 cells?
t11 <- subset(obj, cells = colnames(obj)[obj$paired_1_1 %in% TRUE])
t11$flag <- ifelse(t11$expanded, "expanded clone", "not expanded")
save_plot(DotPlot(t11, features = c("CCL5", "NKG7", "GZMK", "GZMH", "CD8A",
                                    "CCR7", "SELL", "TCF7"), group.by = "flag") +
            RotatedAxis() + labs(x = NULL, y = NULL,
                                 title = "Expanded vs non-expanded TCR 1+1 cells"),
          "16_dotplot_expanded_vs_not", 7, 3.5)

# Cells of expanded clones that fall outside the two test groups. Clones whose
# cells all sit in the CD4 T cluster are plausible CD4 memory expansions; single
# cells with monocyte/NK/pDC labels are more likely doublets or ambient TCR
# (check the dbl_class column).
outside <- obj@meta.data %>%
  filter(expanded, !test_group %in% c("CD8 memory/effector T", "Naive T")) %>%
  select(cell_type, clone_id, clone_size, nFeature_RNA, nCount_RNA, percent.mt,
         any_of(c("dbl_class", "dbl_score")))
write.csv(outside, file.path(RES_DIR, "expanded_cells_outside_test_groups.csv"))
print(table(outside$cell_type))

# ---- 14. Internal consistency checks ----------------------------------------
stopifnot(
  sum(cs) == sum(obj$paired_1_1, na.rm = TRUE),                        # every 1+1 cell has one clone
  sum(by_type$cells_in_expanded) == sum(obj$expanded),                  # no expanded cell lost
  nrow(outside) + sum(res_main$cells_in_expanded) == sum(obj$expanded), # groups + outside = all
  sum(exp_clones$size) == sum(obj$expanded)
)
message("Consistency checks passed.")

# ---- 15. Save ---------------------------------------------------------------
saveRDS(obj, file.path(RES_DIR, "obj_labeled.rds"))
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, "sessionInfo.txt"))
message("Done.")
