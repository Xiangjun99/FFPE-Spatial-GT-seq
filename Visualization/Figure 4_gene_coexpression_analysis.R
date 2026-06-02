```r

library(hdWGCNA)
library(WGCNA)
library(Seurat)
library(tidyverse)
library(igraph)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(stringr)

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 15)

load("/gpfs/ycga/work/liu_yang/xd97/FFPE_triomics_analysis/2_FFPE_RNA_final/12_Bladdercancer_RNA111/visualization/lncRNA/@260329/FinalV1_pbmc_lncRNA_SeuratV5_V2_260329.rdata")

# DimPlot(scRNA, label = TRUE)

Bra_rna <- pbmc

DefaultAssay(Bra_rna) <- "RNA"
Bra_rna <- NormalizeData(Bra_rna)

#===========================================================================================
# 1. Pre-processing
#===========================================================================================

Bra_rna <- SetupForWGCNA(
  Bra_rna,
  wgcna_name = "Bra",
  gene_select = "fraction",
  fraction = 0.05
)

Bra_rna <- MetacellsByGroups(
  seurat_obj = Bra_rna,
  group.by = c("celltype", "orig.ident"),
  k = 25,
  max_shared = 20,
  reduction = "umap",
  ident.group = "celltype"
)

Bra_rna <- NormalizeMetacells(Bra_rna)

#===========================================================================================
# 2. Co-expression network analysis
#===========================================================================================

Bra_rna <- SetDatExpr(Bra_rna, assay = "RNA", slot = "data")

Bra_rna <- TestSoftPowers(Bra_rna, networkType = "signed")

plot <- PlotSoftPowers(Bra_rna)

wrap_plots(plot, ncol = 2)

Bra_rna <- ConstructNetwork(
  Bra_rna,
  soft_power = 4,
  tom_name = "Bra_Test",
  setDatExpr = FALSE
)

PlotDendrogram(Bra_rna, main = "scRNA hdWGCNA Dendrogram")

TOM <- GetTOM(Bra_rna)

#===========================================================================================
# 3. Compute module eigengenes and connectivity
#===========================================================================================

Bra_rna <- ScaleData(Bra_rna)
Bra_rna <- ModuleEigengenes(Bra_rna)

hMEs <- GetMEs(Bra_rna)

MEs <- GetMEs(Bra_rna, harmonized = FALSE)

Bra_rna <- ModuleConnectivity(Bra_rna)

Bra_rna <- ResetModuleNames(Bra_rna, new_name = "M")

modules <- GetModules(Bra_rna)

mod_color_df <- GetModules(Bra_rna) %>%
  dplyr::select(c(module, color)) %>%
  distinct() %>%
  arrange(module)

n_mods <- nrow(mod_color_df) - 1

newcolor <- c(
  "#f4c40f", "#fe9b00", "#d8443c", "#de597c", "#e87b89",
  "#633372", "#1f6e9c", "#2b9b81", "#92c051"
)

Bra_rna <- ResetModuleColors(Bra_rna, newcolor)

modules <- GetModules(Bra_rna)
write.csv(modules, file = "modules.csv")

PlotDendrogram(Bra_rna, main = "scRNA hdWGCNA Dendrogram")

hub_df <- GetHubGenes(Bra_rna, n_hubs = 25)
write.csv(hub_df, file = "hub_df.csv")

library(UCell)

Bra_rna <- ModuleExprScore(
  Bra_rna,
  n_genes = 25,
  method = "UCell"
)

save(Bra_rna, file = "Bra_rna_wgcna.RData")

#===========================================================================================
# 4. Visualize the network
#===========================================================================================

PlotKMEs(Bra_rna, ncol = 3)

plot_hMEs <- ModuleFeaturePlot(
  Bra_rna,
  reduction = "umap",
  features = "hMEs",
  order = TRUE,
  raster = TRUE
)

pdf("5_ModuleFeaturePlot_hMEs.pdf", width = 10, height = 8)
wrap_plots(plot_hMEs, ncol = 3)
dev.off()

plot_score <- ModuleFeaturePlot(
  Bra_rna,
  reduction = "umap",
  features = "hMEs",
  order = TRUE,
  raster = TRUE,
  ucell = TRUE
)

pdf("5_ModuleFeaturePlot_score.pdf", width = 10, height = 8)
wrap_plots(plot_score, ncol = 3)
dev.off()

Bra_rna@meta.data <- cbind(
  Bra_rna@meta.data,
  GetMEs(Bra_rna, harmonized = TRUE)
)

MEs <- GetMEs(Bra_rna, harmonized = TRUE)
mods <- colnames(MEs)
mods <- mods[mods != "grey"]

DotPlot(Bra_rna, features = mods, group.by = "metacell_grouping") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p1 <- VlnPlot(
  Bra_rna,
  features = c("M2"),
  group.by = "celltype",
  pt.size = 0
)

p1 <- p1 + geom_boxplot(width = 0.25, fill = "white")

p1 <- p1 + xlab("") + ylab("hME") + NoLegend()

p1

ModuleNetworkPlot(Bra_rna, outdir = "./Bra_ModuleNetworks_final")

pdf("8-HubGeneNetworkPlot.pdf", width = 8, height = 5)

HubGeneNetworkPlot(
  Bra_rna,
  n_hubs = 10,
  n_other = 5,
  edge_prop = 0.75,
  mods = "all"
)

dev.off()

Bra_rna <- RunModuleUMAP(
  Bra_rna,
  n_hubs = 10,
  n_neighbors = 10,
  min_dist = 0.1
)

pdf("9_ModuleUMAPPlot_fianl.pdf", width = 8, height = 8)

ModuleUMAPPlot(
  Bra_rna,
  edge.alpha = 0.25,
  sample_edges = TRUE,
  edge_prop = 0.1,
  label_hubs = 3,
  keep_grey_edges = FALSE
)

dev.off()

#===========================================================================================
# 5. Enrichment analysis
#===========================================================================================

library(enrichR)

dbs <- c("GO_Biological_Process_2025")

Bra_rna <- RunEnrichr(
  Bra_rna,
  dbs = dbs,
  max_genes = 100
)

enrich_df <- GetEnrichrTable(Bra_rna)
write.csv(enrich_df, "05_enrich_GO_2023.csv")

pdf("10_module_GO.pdf", width = 8, height = 8)

EnrichrDotPlot(
  Bra_rna,
  mods = "all",
  database = "GO_Biological_Process_2025",
  n_terms = 5
)

dev.off()
```
