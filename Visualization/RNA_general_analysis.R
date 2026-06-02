library(Seurat)
library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(dplyr)
library(Matrix)
library(BuenColors)


dir <- "/path/to/spatial-GT-seq/"  
setwd(dir) 


##read expression matrix
my_data <- read.table(file = 'RNA_Filtered_matrix.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
data_filtered <- my_data

data_filtered$X.1 = NULL
count <- rowSums(data_filtered[,2:ncol(data_filtered)], na.rm = TRUE)
data_filtered_binary <- data_filtered[,2:ncol(data_filtered)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)

#extract the coordinates of each pixel
test <- data_filtered %>% separate(X, c("A", "B"),  sep = "x")

df1 <- data.frame(count = count, type = "UMI")
df2 <- data.frame(count = gene_count, type = "Gene")
df_all <- rbind(df1, df2)

x_max <- 20000  

# plot
pdf(file = "UMI_Gene_distribution.pdf", width=8.6, height=8.6)
ggplot(df_all, aes(x = count, fill = type, color = type)) + 
  geom_histogram(aes(y = ..density..), binwidth = x_max / 20, alpha = 0.3, position = "identity") +
  geom_density(alpha = 0.3, size = 1) +
  scale_x_continuous(name = "Count", limits = c(0, x_max)) +
  scale_y_continuous(name = "Density", expand = c(0, 0)) +
  scale_fill_manual(values = c("UMI" = "#1F5188", "Gene" = "#BA383F")) +
  scale_color_manual(values = c("UMI" = "#1F5188", "Gene" = "#BA383F")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.text = element_text(colour = "black", size = 20),
    axis.title = element_text(colour = "black", size = 25, face = "bold"),
    legend.text = element_text(colour = "black", size = 20),
    legend.title = element_text(colour = "black", size = 20, face = "bold"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
    axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid')
  ) +
  labs(title = "UMI and Gene Count Distribution", fill = "Type", color = "Type")
dev.off()

#UMI heatmap
pdf(file = paste("UMI_heatmap.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count))  +
  scale_color_gradientn(colours = jdb_palette("brewer_spectra"),limits=c(0,10000),
                        oob = scales::squish) +
  ggtitle("UMI") +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 16, size = 4)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

#Gene heatmap
pdf(file = paste("Gene_heatmap.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=gene_count))  +
  scale_color_gradientn(colours = jdb_palette("brewer_yes"),limits=c(0,5000),
                        oob = scales::squish) +
  ggtitle("Gene") +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 16, size = 4)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

data1 <- read.table("RNA_Filtered_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)

data3 <- data.frame(X = rownames(data1), data1)
temp1 <- data3 %>% separate(X, c("A", "B"),  sep = "x")
temp1$X.1 = NULL

temp1$A = NULL
temp1$B = NULL
temp1$unmapped = NULL
data2 <- t(temp1)
sample1.name <- "atrium"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)


alphabet <- c('0' = "#73A1BE",
              '1' = "#EB545C",
              '2' = "#efd510",
              '3' = '#7F3C8D',
              '4' = '#80FF08')


pbmc          <- CreateSeuratObject(matrix1.data, min.cells = 1,min.features = 10, project = sample1.name)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
ElbowPlot(pbmc, ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution=0.4, verbose = FALSE)

p1 <- DimPlot(pbmc, label = TRUE) + NoLegend() +
  scale_color_manual(values = alphabet[1:(pbmc@active.ident %>% unique %>% length )])
p1
ggsave("UMAP.pdf", plot = p1, path = dir, dpi = 600, bg = NULL, height = 8, width = 8)

ident <- Idents(pbmc)
df <- data.frame(ident[])
df1 <-data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

pdf(file = paste("Spatial_cluster.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) + 
  scale_color_manual(values = alphabet[1:(pbmc@active.ident %>% unique %>% length )]) + 
  ggtitle("UMAP") +
  geom_point(shape = 16, size = 4.5)+
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

levels(pbmc)
pbmc.markers<- FindAllMarkers(pbmc,logfc.threshold = 0)
pbmc <-ScaleData(object=pbmc, features = rownames(pbmc))
pbmc.markers.top5<- pbmc.markers %>% group_by(cluster) %>% top_n(wt = avg_log2FC , n =  10)
plotE <- DoHeatmap(pbmc, features = pbmc.markers.top5$gene, group.by = "seurat_clusters",  
                   group.colors = c('0' = "#73A1BE",
                                    '1' = "#EB545C",
                                    '2' = "#efd510",
                                    '3' = '#7F3C8D',
                                    '4' = '#80FF08'))+ 
  theme(axis.text=element_text(size=5))+ 
  scale_fill_gradientn(colors=c("#192E5B","#1D65A6","#72A2C0","#F3E96B","#F2A104"))
plotE

saveRDS(pbmc.markers,file="DEGs.rds")
write.table(pbmc.markers,file='DEGs.txt')
ggsave("HEATMAP_top_markergene.pdf", plot = plotE, path = dir, dpi = 600, bg = NULL, height = 6, width = 8)
save(pbmc, file= "pbmc.rdata")
