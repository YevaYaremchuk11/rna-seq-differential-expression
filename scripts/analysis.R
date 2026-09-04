# RNA-seq differential expression analysis of GSE213111
# Endothelial inflammatory response
#
# Dataset: GSE213111
# Organism: Homo sapiens
# Main tools: GEOquery, DESeq2, pheatmap, AnnotationDbi,
#             org.Hs.eg.db, clusterProfiler, enrichplot, ggplot2
#
# DEG criteria used in this analysis:
#   adjusted p-value < 0.05
#   absolute log2 fold change >= 1

# -----------------------------
# 1. Packages and output folders
# -----------------------------

library(BiocManager)
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Download GEO metadata
# -----------------------------

gse <- getGEO("GSE213111", GSEMatrix = FALSE)

# Extract sample-level treatment information and descriptions.
treatments <- unlist(lapply(
  gse@gsms,
  function(x) x@header$characteristics_ch1[[4]]
))

desc <- unlist(lapply(
  gse@gsms,
  function(x) x@header$description
))

# -----------------------------
# 3. Download count matrix
# -----------------------------

counts_file <- "GSE213111_ARH01_counts.txt.gz"

download.file(
  gse@header$supplementary_file,
  destfile = counts_file
)

header <- strsplit(
  readLines(counts_file, n = 2)[2],
  "\t"
)[[1]]

sample_names <- basename(header[7:84])
sample_names <- gsub(
  "_Aligned.sortedByCoord.out.bam",
  "",
  sample_names
)

counts <- read.table(
  counts_file,
  sep = "\t",
  header = TRUE
)

counts_matrix <- counts[, c(1, 7:84)]
colnames(counts_matrix)[2:79] <- sample_names

# -----------------------------
# 4. Build sample metadata
# -----------------------------

gsm_ordered <- names(desc)[match(sample_names, desc)]

treatments_ordered <- treatments[
  match(gsm_ordered, names(treatments))
]

metadata <- data.frame(
  sample = sample_names,
  gsm = gsm_ordered,
  treatments = treatments_ordered
)

metadata$treatment_clean <- gsub(
  "treatment: ",
  "",
  metadata$treatments
)

metadata$time <- sub(
  ".* ",
  "",
  metadata$treatment_clean
)

metadata$treatment <- sub(
  " .*",
  "",
  metadata$treatment_clean
)

metadata$treatments <- NULL
metadata$treatment_clean <- NULL

# SILAC labels are balanced across stimulated treatment/time groups
# and are treated here as replicate labels rather than a main design factor.
metadata$silac <- sub(
  ".*([HLM]).*",
  "\\1",
  metadata$sample
)

# Convert categorical variables to factors.
metadata$treatment <- factor(metadata$treatment)
metadata$time <- factor(metadata$time)
metadata$silac <- factor(metadata$silac)

metadata$treatment <- relevel(
  metadata$treatment,
  ref = "Sham"
)

# Combine treatment and time into one factor.
metadata$group <- paste(
  metadata$treatment,
  metadata$time,
  sep = "_"
)

metadata$group <- factor(metadata$group)
metadata$group <- relevel(
  metadata$group,
  ref = "Sham_0min"
)

# Confirm count matrix and metadata are aligned.
stopifnot(identical(
  colnames(counts_matrix)[2:79],
  metadata$sample
))

# -----------------------------
# 5. Filter low-count genes
# -----------------------------

# Keep genes with count >= 10 in at least 3 samples.
gene_n_samples <- rowSums(counts_matrix[, 2:79] >= 10)
keep <- gene_n_samples >= 3

counts_filtered <- counts_matrix[keep, ]

stopifnot(anyDuplicated(counts_filtered$Geneid) == 0)

rownames(counts_filtered) <- counts_filtered$Geneid
counts_filtered$Geneid <- NULL

counts_matrix_deseq <- as.matrix(counts_filtered)

# -----------------------------
# 6. DESeq2 analysis
# -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix_deseq,
  colData = metadata,
  design = ~ group
)

dds <- DESeq(dds)

# Available contrasts.
contrasts <- resultsNames(dds)[-1]

# -----------------------------
# 7. Differential expression for all contrasts
# -----------------------------

deg_results <- list()

for (name in contrasts) {
  res <- results(dds, name = name)

  deg <- res[
    !is.na(res$padj) &
      res$padj < 0.05 &
      abs(res$log2FoldChange) >= 1,
  ]

  deg_results[[name]] <- deg
}

# Summary of DEG counts.
deg_summary <- data.frame(
  contrast = names(deg_results),
  n_DEGs = sapply(deg_results, nrow)
)

deg_up_down <- lapply(deg_results, function(res) {
  c(
    up = sum(res$log2FoldChange >= 1, na.rm = TRUE),
    down = sum(res$log2FoldChange <= -1, na.rm = TRUE)
  )
})

deg_summary$up <- sapply(deg_up_down, function(x) x["up"])
deg_summary$down <- sapply(deg_up_down, function(x) x["down"])

deg_summary$time <- sub(
  ".*group_.*_(2min|5min|10min|30min|4h|8h|12h|24h)_vs_.*",
  "\\1",
  deg_summary$contrast
)

deg_summary$treatment <- sub(
  "group_(TNF\\.IFNgamma|IFNgamma|TNF)_.*",
  "\\1",
  deg_summary$contrast
)

deg_summary$treatment[
  deg_summary$treatment == "TNF.IFNgamma"
] <- "TNF+IFNgamma"

deg_summary$time <- factor(
  deg_summary$time,
  levels = c(
    "2min", "5min", "10min", "30min",
    "4h", "8h", "12h", "24h"
  )
)

deg_plot <- deg_summary[, c(
  "treatment", "time", "n_DEGs", "up", "down"
)]

write.csv(
  deg_summary,
  file.path("results/tables", "DEG_summary.csv"),
  row.names = FALSE
)

# Store complete DESeq2 result tables.
all_results <- lapply(contrasts, function(name) {
  results(dds, name = name)
})

names(all_results) <- contrasts

# -----------------------------
# 8. DEG trajectory over time
# -----------------------------

deg_trajectory_plot <- ggplot(
  deg_plot,
  aes(
    x = time,
    y = n_DEGs,
    group = treatment,
    color = treatment
  )
) +
  geom_line() +
  geom_point() +
  labs(
    x = "Time",
    y = "Number of DEGs",
    title = "Differentially expressed genes over time"
  ) +
  theme_minimal()

ggsave(
  "figures/DEG_over_time.png",
  deg_trajectory_plot,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------
# 9. TNF 4h volcano plot
# -----------------------------

res_TNF_4h <- all_results[[
  "group_TNF_4h_vs_Sham_0min"
]]

plot_data <- as.data.frame(res_TNF_4h)
plot_data$gene <- rownames(plot_data)
plot_data$significance <- "Not significant"

plot_data$significance[
  !is.na(plot_data$padj) &
    plot_data$padj < 0.05 &
    abs(plot_data$log2FoldChange) >= 1
] <- "DEG"

volcano_plot <- ggplot(
  plot_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(alpha = 0.5) +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    title = "TNF 4h vs Sham 0min"
  ) +
  theme_minimal()

ggsave(
  "figures/volcano_TNF_4h.png",
  volcano_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 10. PCA and heatmap QC
# -----------------------------

vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(
  vsd,
  intgroup = c("treatment", "time")
)

ggsave(
  "figures/PCA.png",
  pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

top_genes <- order(
  apply(assay(vsd), 1, var),
  decreasing = TRUE
)[1:min(50, nrow(assay(vsd)))]

heatmap_data <- assay(vsd)[top_genes, ]

annotation_col <- metadata[, c("treatment", "time")]
rownames(annotation_col) <- metadata$sample

pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Top 50 most variable genes",
  filename = "figures/heat.png",
  width = 8,
  height = 8
)

# -----------------------------
# 11. Gene annotation for TNF 4h
# -----------------------------

deg_TNF_4h <- res_TNF_4h[
  !is.na(res_TNF_4h$padj) &
    res_TNF_4h$padj < 0.05 &
    abs(res_TNF_4h$log2FoldChange) >= 1,
]

deg_TNF_4h_annotated <- as.data.frame(deg_TNF_4h)
deg_TNF_4h_annotated$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(deg_TNF_4h_annotated),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

deg_TNF_4h_annotated <- deg_TNF_4h_annotated[
  order(deg_TNF_4h_annotated$padj),
]

write.csv(
  deg_TNF_4h_annotated,
  file.path("results/tables", "DEGs_TNF_4h_vs_Sham_0min_annotated.csv")
)

# -----------------------------
# 12. GO enrichment at 4h
# -----------------------------

get_gene_symbols <- function(deg_table) {
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(deg_table),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )

  symbols[!is.na(symbols)]
}

# TNF 4h
genes_TNF_4h <- get_gene_symbols(
  deg_results[["group_TNF_4h_vs_Sham_0min"]]
)

ego_TNF_4h <- enrichGO(
  gene = genes_TNF_4h,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# IFNgamma 4h
genes_IFNgamma_4h <- get_gene_symbols(
  deg_results[["group_IFNgamma_4h_vs_Sham_0min"]]
)

ego_IFNgamma_4h <- enrichGO(
  gene = genes_IFNgamma_4h,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# TNF + IFNgamma 4h
genes_TNF_IFNgamma_4h <- get_gene_symbols(
  deg_results[["group_TNF.IFNgamma_4h_vs_Sham_0min"]]
)
ego_TNF_IFNgamma_4h <- enrichGO(
  gene = genes_TNF_IFNgamma_4h,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Save GO tables.
write.csv(
  as.data.frame(ego_TNF_4h),
  file.path("results/tables", "GO_TNF_4h.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_IFNgamma_4h),
  file.path("results/tables", "GO_IFNgamma_4h.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_TNF_IFNgamma_4h),
  file.path("results/tables", "GO_TNF_IFNgamma_4h.csv"),
  row.names = FALSE
)

# Save GO dot plots.
ggsave(
  "figures/GO_TNF_4h.png",
  dotplot(ego_TNF_4h, showCategory = 15) + ggtitle("GO enrichment: TNF 4h"),
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GO_IFNgamma_4h.png",
  dotplot(ego_IFNgamma_4h, showCategory = 15) + ggtitle("GO enrichment: IFNgamma 4h"),
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GO_TNF_IFNgamma_4h.png",
  dotplot(ego_TNF_IFNgamma_4h, showCategory = 15) + ggtitle("GO enrichment: TNF + IFNgamma 4h"),
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# 13. Shared GO processes at 4h
# -----------------------------

common_GO_4h <- Reduce(
  intersect,
  list(
    ego_TNF_4h@result$Description,
    ego_IFNgamma_4h@result$Description,
    ego_TNF_IFNgamma_4h@result$Description
  )
)

# Merge enrichment statistics for terms shared by all three conditions.
go_common_4h <- merge(
  ego_TNF_4h@result[, c(
    "Description", "FoldEnrichment", "p.adjust"
  )],
  ego_IFNgamma_4h@result[, c(
    "Description", "FoldEnrichment", "p.adjust"
  )],
  by = "Description",
  suffixes = c("_TNF", "_IFNgamma")
)

go_common_4h <- merge(
  go_common_4h,
  ego_TNF_IFNgamma_4h@result[, c(
    "Description", "FoldEnrichment", "p.adjust"
  )],
  by = "Description"
)

go_common_4h <- go_common_4h[
  order(go_common_4h$p.adjust),
]

write.csv(
  go_common_4h,
  file.path("results/tables", "GO_shared_4h.csv"),
  row.names = FALSE
)

# Heatmap of the 20 strongest shared GO processes.
if (nrow(go_common_4h) > 0) {
  go_top20 <- go_common_4h[
    seq_len(min(20, nrow(go_common_4h))),
  ]

  heatmap_go <- go_top20[, c(
    "FoldEnrichment_TNF",
    "FoldEnrichment_IFNgamma",
    "FoldEnrichment"
  )]

  colnames(heatmap_go) <- c(
    "TNF",
    "IFNgamma",
    "TNF+IFNgamma"
  )

  rownames(heatmap_go) <- go_top20$Description

  pheatmap(
    heatmap_go,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Shared GO enrichment at 4h",
    filename = "figures/GO_shared_4h.png",
    width = 9,
    height = 8
  )
}

# -----------------------------
# 14. Top genes at 4h
# -----------------------------

top_genes_4h <- lapply(
  c(
    "group_TNF_4h_vs_Sham_0min",
    "group_IFNgamma_4h_vs_Sham_0min",
    "group_TNF.IFNgamma_4h_vs_Sham_0min"
  ),
  function(x) {
    res <- all_results[[x]]

    res <- res[
      !is.na(res$padj) &
        res$padj < 0.05 &
        abs(res$log2FoldChange) >= 1,
    ]

    res <- res[order(res$padj), ]
    res[seq_len(min(20, nrow(res))), , drop = FALSE]
  }
)

names(top_genes_4h) <- c(
  "TNF",
  "IFNgamma",
  "TNF_IFNgamma"
)

top_genes_4h_symbols <- lapply(
  top_genes_4h,
  function(x) {
    x$gene_symbol <- mapIds(
      org.Hs.eg.db,
      keys = rownames(x),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    x
  }
)

# Genes shared by the top 20 in all three conditions.
common_top_genes_4h <- Reduce(
  intersect,
  lapply(
    top_genes_4h_symbols,
    function(x) unique(na.omit(x$gene_symbol))
  )
)

write.csv(
  data.frame(gene_symbol = common_top_genes_4h),
  file.path("results/tables", "common_top20_genes_4h.csv"),
  row.names = FALSE
)

# Save session information for reproducibility.
writeLines(
  capture.output(sessionInfo()),
  file.path("results/tables", "sessionInfo.txt")
)

# End of analysis.
