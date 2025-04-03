# Source from:
# https://bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/analysis.html

## ----echo = FALSE, message = FALSE, warning = FALSE---------------------------
library(patchwork)
library(cowplot)

## ----load-libs, message = FALSE,  warning = FALSE-----------------------------
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)

## ----eh, message = FALSE------------------------------------------------------
library(ExperimentHub)
eh <- ExperimentHub()

## ----load-data, message = FALSE-----------------------------------------------
sce <- eh[["EH2259"]]
dim(sce)

## ----fil----------------------------------------------------------------------
# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

## ----qc, message = FALSE------------------------------------------------------
# calculate per-cell quality control (QC) metrics
library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

## ----norm---------------------------------------------------------------------
# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
dim(sce)
## ----vst, eval = FALSE--------------------------------------------------------
#  library(sctransform)
#  assays(sce)$vstresiduals <- vst(counts(sce), verbosity = FALSE)$y

## ----prep-sce-----------------------------------------------------------------
sce$id <- paste0(sce$stim, sce$ind)
sce <- prepSCE(sce,
                kid = "cell", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE)  # drop all other colData columns

## ----ids----------------------------------------------------------------------
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

## ----ncells, size = "small"---------------------------------------------------
# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

## ----umap---------------------------------------------------------------------
# compute UMAP using 1st 20 PCs
#sce <- runUMAP(sce, pca = 20)

## -----------------------------------------------------------------------------
# wrapper to prettify reduced dimension plots
# .plot_dr <- function(sce, dr, col)
#   plotReducedDim(sce, dimred = dr, colour_by = col) +
#   guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   theme_minimal() + theme(aspect.ratio = 1)

## ----eval = FALSE-------------------------------------------------------------
#  # downsample to max. 100 cells per cluster
#  cs_by_k <- split(colnames(sce), sce$cluster_id)
#  cs100 <- unlist(sapply(cs_by_k, function(u)
#    sample(u, min(length(u), 100))))
#
#  # plot t-SNE & UMAP colored by cluster & group ID
#  for (dr in c("TSNE", "UMAP"))
#    for (col in c("cluster_id", "group_id"))
#      .plot_dr(sce[, cs100], dr, col)

## ----dr-ids, echo = FALSE, results = "asis", fig.height = 4, fig.width = 12, fig.cap = "Dimension reduction plots. Cells are colored by cluster ID (A) and group ID (B), respectively. For each cluster, at most 100 cells were sampled for plotting."----
# cs_by_k <- split(colnames(sce), sce$cluster_id)
# cs100 <- unlist(sapply(cs_by_k, function(u)
#   sample(u, min(length(u), 100))))

# for (dr in c("TSNE", "UMAP")) {
#   cat("#### ", dr, "{-}\n")
#   ps <- lapply(c("cluster_id", "group_id"),
#                function(col) .plot_dr(sce[, cs100], dr, col = col))
#   assign(paste0("ps_", tolower(dr)), ps)
#   print(wrap_plots(ps, align = "vh", labels = c("A", "B")))
#   cat("\n\n")
# }

## ----echo = FALSE, out.height = 4, fig.cap = "Schematic overview of cell- and sample-level approaches for DS analysis. Top panels show a schematic of the data distributions or aggregates across samples (each violin is a group or sample; each dot is a sample) and conditions (blue or orange). The bottom panels highlight the data organization in sub-matrix slices of the original count table."----
#knitr::include_graphics(system.file('extdata', '1d.png', package = 'muscat'))

## ----agg----------------------------------------------------------------------
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)
# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## ----pb-mds, fig.height = 4, fig.cap = "Pseudobulk-level multidimensional scaling (MDS) plot. Each point represents a cluster-sample instance; points are colored by cluster ID and shaped by group ID."----
pb_mds <- pbMDS(pb)

## ----message = FALSE, fig.height = 4, fig.cap = "Pseudobulk-level MDS plot v2. Default plotting aesthetics were modified to change shaping of groups, coloring of clusters, as well as point size and transparency."----
# use very distinctive shaping of groups & change cluster colors
# pb_mds <- pb_mds +
#   scale_shape_manual(values = c(17, 4)) +
#   scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
# # change point size & alpha
# pb_mds$layers[[1]]$aes_params$size <- 5
# pb_mds$layers[[1]]$aes_params$alpha <- 0.6
# pb_mds

## -----------------------------------------------------------------------------
# run DS analysis
res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)
# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))

## ----eval = FALSE-------------------------------------------------------------
#  # construct design & contrast matrix
#  ei <- metadata(sce)$experiment_info
#  mm <- model.matrix(~ 0 + ei$group_id)
#  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
#  contrast <- makeContrasts("stim-ctrl", levels = mm)
#
#  # run DS analysis
#  pbDS(pb, design = mm, contrast = contrast)

## ----mm, eval = FALSE---------------------------------------------------------
#  # 1st approach
#  mm <- mmDS(sce, method = "dream",
#    n_cells = 10, n_samples = 2,
#    min_counts = 1, min_cells = 20)
#
#  # 2nd & 3rd approach
#  mm <- mmDS(sce, method = "vst", vst = "sctransform")
#  mm <- mmDS(sce, method = "nbinom")

## -----------------------------------------------------------------------------
# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# view top 2 hits in each cluster
top2 <- bind_rows(lapply(tbl_fil, top_n, 2, p_adj.loc))
format(top2[, -ncol(top2)], digits = 2)

## ----frq----------------------------------------------------------------------
frq <- calcExprFreqs(sce, assay = "counts", th = 0)
# one sheet per cluster
assayNames(frq)
# expression frequencies in each
# sample & group; 1st cluster
t(head(assay(frq), 5))

## -----------------------------------------------------------------------------
gids <- levels(sce$group_id)
frq10 <- vapply(as.list(assays(frq)),
                function(u) apply(u[, gids] > 0.1, 1, any),
                logical(nrow(sce)))
t(head(frq10))

tbl_fil2 <- lapply(kids, function(k)
  dplyr::filter(tbl_fil[[k]],
                gene %in% names(which(frq10[, k]))))

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil2, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  # tidy format; attach pre-computed expression frequencies
#  resDS(sce, res, bind = "row", frq = frq)
#
#  # big-table (wide) format; attach CPMs
#  resDS(sce, res, bind = "col", cpm = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  # compute expression frequencies on the fly
#  resDS(sce, res, frq = TRUE)

## ----upset, fig.width = 10, fig.cap = "Upset plot. Included are DS findings (FDR < 0.05, |logFC| > 1) across all clusters; shown are the 50 most frequent interactions."----
library(UpSetR)
de_gs_by_k <- map(tbl_fil, "gene")
upset(fromList(de_gs_by_k))

## ----fig.width = 14, fig.height = 8, fig.cap = "t-SNE colored by gene expression. Show are t-SNE projections with cells colored by the expression of the top-8 DS genes. For each cluster, at most 100 cells were sampled for plotting."----
# pull top-8 DS genes across all clusters
top8 <- bind_rows(tbl_fil) %>%
  slice_min(p_adj.loc, n = 8,
            with_ties = FALSE) %>%
  pull("gene")

# for ea. gene in 'top8', plot t-SNE colored by its expression
# ps <- lapply(top8, function(g)
#   .plot_dr(sce[, cs100], "TSNE", g) +
#     ggtitle(g) + theme(legend.position = "none"))

# arrange plots
# plot_grid(plotlist = ps, ncol = 4, align = "vh")

## ----violins, fig.width = 10, fig.height = 5, fig.cap = "Violin plots. Show are the top 6 hits (lowest adj. p-value) for the B cells cluster. Each violin is a sample; points are colored by group ID."----
plotExpression(sce[, sce$cluster_id == "B cells"],
               features = tbl_fil$`B cells`$gene[seq_len(6)],
               x = "sample_id", colour_by = "group_id", ncol = 3) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----pb-hm-1------------------------------------------------------------------
# top-5 DS genes per cluster
pbHeatmap(sce, res, top_n = 5)

## ----pb-hm-2------------------------------------------------------------------
# top-20 DS genes for single cluster
pbHeatmap(sce, res, k = "B cells")

## ----pb-hm-3------------------------------------------------------------------
# single gene across all clusters
pbHeatmap(sce, res, g = "ISG20")

## ----session-info-------------------------------------------------------------
sessionInfo()
