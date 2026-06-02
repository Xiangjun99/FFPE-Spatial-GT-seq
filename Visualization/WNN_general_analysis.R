library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)

# -----------------------------
# 1. load
# -----------------------------
read_cell_by_feature <- function(file) {
  dt <- fread(file)
  
  cell_ids <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1, with = FALSE])
  
  rownames(mat) <- cell_ids
  colnames(mat) <- make.unique(colnames(mat))
  storage.mode(mat) <- "numeric"
  
  mat <- t(mat)
  mat <- Matrix(mat, sparse = TRUE)
  
  return(mat)
}

rna <- read_cell_by_feature("RNA_Filtered_matrix_correct.tsv")
adt <- read_cell_by_feature("ADT_Filtered_matrix_correct.tsv")

cat("RNA matrix:", nrow(rna), "features x", ncol(rna), "cells\n")
cat("ADT matrix:", nrow(adt), "features x", ncol(adt), "cells\n")

# -----------------------------
# 2.find common cells
# -----------------------------
common_cells <- intersect(colnames(rna), colnames(adt))

cat("Common cells:", length(common_cells), "\n")

if (length(common_cells) == 0) {
  stop("RNA and ADT no cell ID.")
}

rna <- rna[, common_cells, drop = FALSE]
adt <- adt[, common_cells, drop = FALSE]

# -----------------------------
# 3. build Seurat object
# -----------------------------
obj <- CreateSeuratObject(
  counts = rna,
  assay = "RNA",
  project = "RNA_ADT_WNN",
  min.cells = 1,
  min.features = 10
)

adt <- adt[, Cells(obj), drop = FALSE]

obj[["ADT"]] <- CreateAssayObject(counts = adt)

# -----------------------------
# 4. RNA preprocessing
# -----------------------------

DefaultAssay(obj) <- "RNA"

obj <- SCTransform(
  obj,
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = FALSE
)

obj <- RunPCA(
  obj,
  assay = "SCT",
  reduction.name = "rna.sct.pca",
  reduction.key = "rnaSCTPC_",
  npcs = 50,
  verbose = FALSE
)

# -----------------------------
# 5. ADT preprocessing
# -----------------------------
adt_features <- rownames(obj[["ADT"]])

obj <- SCTransform(
  obj,
  assay = "ADT",
  new.assay.name = "ADT_SCT",
  residual.features = adt_features,
  variable.features.n = length(adt_features),
  verbose = FALSE
)

DefaultAssay(obj) <- "ADT_SCT"
VariableFeatures(obj) <- rownames(obj[["ADT_SCT"]])

n_adt_pcs <- min(
  18,
  length(VariableFeatures(obj)) - 1,
  ncol(obj) - 1
)

if (n_adt_pcs < 2) {
  stop("ADT features too less, cannot do PCA/WNN。")
}

obj <- RunPCA(
  obj,
  assay = "ADT_SCT",
  features = VariableFeatures(obj),
  reduction.name = "adt.sct.pca",
  reduction.key = "adtSCTPC_",
  npcs = n_adt_pcs,
  verbose = FALSE
)


# -----------------------------
# 6. WNN analysis
# -----------------------------
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("rna.sct.pca", "adt.sct.pca"),
  dims.list = list(1:10, 1:5),
  modality.weight.name = "RNA.weight"
)

obj <- RunUMAP(
  obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

obj <- FindClusters(
  obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 1.0
)
# -----------------------------
# 7. visualization
# -----------------------------
DimPlot(obj, reduction = "wnn.umap", label = TRUE) + NoLegend()

FeaturePlot(
  obj,
  reduction = "wnn.umap",
  features = "RNA.weight"
)

# save
save(obj, file="seurat_WNN_RNA-SCT_ADT-SCT.rdata")


###check modality in each cluster
meta <- obj@meta.data


meta$cluster <- as.factor(Idents(obj))


meta$RNA.weight <- obj$RNA.weight


meta$ADT.weight <- 1 - meta$RNA.weight

cluster_weight_summary <- meta %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    RNA_mean = mean(RNA.weight, na.rm = TRUE),
    ADT_mean = mean(ADT.weight, na.rm = TRUE),
    RNA_median = median(RNA.weight, na.rm = TRUE),
    ADT_median = median(ADT.weight, na.rm = TRUE),
    .groups = "drop"
  )

cluster_weight_summary

cluster_weight_long <- cluster_weight_summary %>%
  dplyr::select(cluster, RNA_mean, ADT_mean) %>%
  pivot_longer(
    cols = c(RNA_mean, ADT_mean),
    names_to = "modality",
    values_to = "weight"
  ) %>%
  mutate(
    modality = recode(
      modality,
      RNA_mean = "RNA",
      ADT_mean = "ADT"
    )
  )

p_weight_bar <- ggplot(cluster_weight_long, aes(x = cluster, y = weight, fill = modality)) +
  geom_col(width = 0.8) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Mean modality weight",
    fill = "Modality"
  ) +
  ylim(0, 1)

p_weight_bar


# -----------------------------
# 8. Cluster Spatial Plot
# -----------------------------

library(ggplot2)

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

coords <- do.call(rbind, strsplit(meta$cell_id, "x"))
meta$spatial_col <- as.numeric(coords[, 1])   # A
meta$spatial_row <- as.numeric(coords[, 2])   # B

head(meta[, c("cell_id", "spatial_col", "spatial_row", "seurat_clusters")])

meta$seurat_clusters <- as.factor(meta$seurat_clusters)


#nclust <- length(unique(meta$seurat_clusters))
alphabet <- c('0' = "#63B5B7",  # aqua green
              '1' = "#F9D367",  # warm gold
              '2' = "#C2412D",  # deep brick red
              '3' = "#FB9A99",  # muted blue
              '4' = "#8C6BB1",  # dusty purple
              '5' = "#4C6A9A",  # olive green
              '6' = "#FF9900",  # dusty purple
              '7' = "#85C17E",
              
              '8' = "#289E92",  # muted blue
              '9' = "#DCCD58",  # dusty purple
              '10' = "#CB7E83",  # olive green
              '11' = "#C24976",  # dusty purple
              '12' = "#FC6FCF"
)

cols <- alphabet[1:nclust]


p_spatial <- ggplot(meta, aes(x = spatial_col, y = spatial_row, color = seurat_clusters)) +
  geom_point(shape = 16, size = 4.5) +
  scale_color_manual(values = cols) +
  ggtitle("Spatial clusters") +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(
    name = "X",
    limits = c(NA, NA),
    expand = expansion(mult = c(-0.013, -0.013))
  ) +
  scale_y_reverse(
    name = "Y",
    limits = c(NA, NA),
    expand = expansion(mult = c(-0.013, 0.008))
  ) +
  coord_equal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold", color = "white"),
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "black", color = "black")
  )

p_spatial

ggsave("WNN_cluster_spatial.pdf", p_spatial, width = 11.6, height = 11.6)


# ------------------------------------
# 9. RNA gene expression spatial plot
# ------------------------------------
gene <- "CD3D"

plot_df <- FetchData(obj, vars = gene)
plot_df$cell_id <- rownames(plot_df)

coords <- do.call(rbind, strsplit(plot_df$cell_id, "x"))
plot_df$spatial_col <- as.numeric(coords[, 1])
plot_df$spatial_row <- as.numeric(coords[, 2])

ggplot(plot_df, aes(x = spatial_col, y = spatial_row, fill = .data[[gene]])) +
  geom_tile() +
  coord_fixed() +
  scale_y_reverse() +
  theme_classic() +
  labs(fill = gene)


# ------------------------------------
# 10. ADT protein spatial plot
# ------------------------------------

DefaultAssay(obj) <- "ADT"

adt_df <- FetchData(obj, vars = "CD274")
adt_df$cell_id <- rownames(adt_df)

coords <- do.call(rbind, strsplit(adt_df$cell_id, "x"))
adt_df$spatial_col <- as.numeric(coords[, 1])
adt_df$spatial_row <- as.numeric(coords[, 2])

ggplot(adt_df, aes(x = spatial_col, y = spatial_row, fill = CD274)) +
  geom_tile() +
  coord_fixed() +
  scale_y_reverse() +
  theme_classic() +
  labs(fill = "CD274")
