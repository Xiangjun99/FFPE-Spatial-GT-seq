setwd("/path/to/spatial-GT-seq/")

library(copykit)
library(SummarizedExperiment)
library(magrittr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(cowplot)
library(GenomicRanges)
library(BiocParallel)

register(MulticoreParam(workers = 32, progressbar = T),default = T)

copykit <-runVarbin('bamFile',
                    genome='hg38', 
                    resolution = "1Mb",
                    remove_Y=T,
                    remove_X = TRUE,
                    is_paired_end=T,
                    min_bincount=1)

save(copykit, file = "copykit_V1.rdata")

col_fun=circlize::colorRamp2(breaks = c(-0.5,0,0.5), 
                             c("#74add1","white","#fd7f64"))

copykit <- runMetrics(copykit)
copykit <- findAneuploidCells(copykit, remove_XY = T, simul = F)

copykit <- runUmap(copykit)

copykit <- findSuggestedK(copykit)
p1 <- plotSuggestedK(copykit)
p1

subclone_pal = c(
  'c1' = "#EB545C",
  'c2' = "#5085C4",
  'c3' = "#F6EC1B",
  'c4' = '#87C55F',
  'c5' = '#CC66FF',
  'c6' = '#59A14F')

copykit  <- findClusters(copykit, k_subclones = 20, method ="leiden") #"hdbscan", "leiden", "louvain"

p2 <- plotUmap(copykit, label = 'subclones') + 
  scale_fill_manual(values = subclone_pal, name = "subclones")
p2
ggsave("CNV_UMAP.pdf", plot = p2, width = 8, height = 8, dpi = 300)

copykit <- calcConsensus(copykit)
copykit <- runConsensusPhylo(copykit)

pdf("CNV_heatmap.pdf", width = 16, height = 8)
plotHeatmap(
  copykit,
  genome = "hg38",
  consensus = TRUE,
  label = c('subclones'),
  order_cells = 'consensus_tree',
  col = col_fun,
  label_colors = list(subclones = subclone_pal),
  n_threads = 100
)
dev.off()

#spatial plot
A <- data.frame(copykit@colData@listData[["sample"]])
colnames(A)[1] <- "cell"
A$ident <- copykit@colData@listData[["subclones"]]
A$sample <- sub("_.*", "", A$cell)

A[, 1] <- gsub("^sort_", "", A[, 1])  
A[, 1] <- gsub("\\+", "", A[, 1])  

#change cell barcode to spatial barcode
lookup = read.table("lookup.csv", as.is = TRUE, sep = ",")
temp = lookup[match(A$cell,lookup$V1),"V2"]
A$name =temp

#add position information to remove the cell outlier of tissue
location <- read.table("tissue_position.txt", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])
x = x[-1]

A_filtered <- A[A$name %in% x,]

test <- A_filtered %>% separate(name, c("A", "B"),  sep = "x")

pdf(file = paste("Spatial_CNV_clustering.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=ident)) + 
  scale_color_manual(values = subclone_pal)  + 
  ggtitle("UMAP") +
  geom_point(shape = 16, size = 4.0)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) +
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +   
  theme(plot.background = element_rect(fill = "black"))
dev.off()

save(copykit, file = "copykit_V2.rdata")
