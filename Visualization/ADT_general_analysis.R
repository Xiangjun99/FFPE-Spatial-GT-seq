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
library(BuenColors)
library(Matrix)

data_dir <- "/gpfs/ycga/work/liu_yang/xd97/FFPE_triomics_analysis/3_FFPE_protein/7_BladderCancer/ADT312_visualization"  
#dir = "C:/data/Sansing/DIS/DISADT"
setwd(data_dir)

mat <- readMM(file.path(data_dir, "umi_count", "matrix.mtx.gz"))

barcodes <- data.table::fread(
  file.path(data_dir, "umi_count", "barcodes.tsv.gz"),
  header = FALSE
)$V1

features_raw <- readLines(
  gzfile(file.path(data_dir, "umi_count", "features.tsv.gz"))
)
feature_names <- sub("-.*$", "", features_raw)
feature_names <- sub(" \\(.*\\)", "", feature_names)

rownames(mat) <- feature_names
colnames(mat) <- barcodes

mat <- as.matrix(mat)
mat <- t(mat)
mat <- as.data.frame(mat)
mat$cell <- rownames(mat)

lookup = read.table("lookup_72x72.csv", as.is = TRUE, sep = ",")
lookup <- lookup[-1, ]
temp = lookup[match(mat$cell,lookup$V1),"V2"]
mat$name =temp

mat$cell <- NULL

write.table(mat, file = 'expression_matrix_name.tsv', sep = '\t',col.names=TRUE, row.names = FALSE,quote = FALSE)

data1  <- read.table(file = 'expression_matrix_name.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)

#extract the coordinates of each pixel
temp1 <- data1 %>% separate(name, c("A", "B"),  sep = "x")

#repair the strips of the data, here is the col = 34
#the repair was done by averaging the two neighboring columns.
col1 = 4
col2 = 45
col3 = 37
col4 = 33


row1 = 4
row2 = 19
#row3 = 44

# #repair Row
for (i in 1:50) {
  temp1[(temp1$A==i&temp1$B==row1),] = as.integer((as.integer(temp1[(temp1$A==i&temp1$B==row1+1),]) + as.integer(temp1[(temp1$A==i&temp1$B==row1-1),])) / 2)
  temp1[(temp1$A==i&temp1$B==row2),] = as.integer((as.integer(temp1[(temp1$A==i&temp1$B==row2+1),]) + as.integer(temp1[(temp1$A==i&temp1$B==row2-1),])) / 2)
  #temp1[(temp1$A==i&temp1$B==row3),] = as.integer((as.integer(temp1[(temp1$A==i&temp1$B==row3+1),]) + as.integer(temp1[(temp1$A==i&temp1$B==row3-1),])) / 2)
}

# #repair column
for (i in 1:50) {
  temp1[(temp1$A==col1&temp1$B==i),] = as.integer((as.integer(temp1[(temp1$A==(col1-1)&temp1$B==i),]) + as.integer(temp1[(temp1$A==(col1+1)&temp1$B==i),])) / 2)
  temp1[(temp1$A==col2&temp1$B==i),] = as.integer((as.integer(temp1[(temp1$A==(col2-1)&temp1$B==i),]) + as.integer(temp1[(temp1$A==(col2+1)&temp1$B==i),])) / 2)
  temp1[(temp1$A==col3&temp1$B==i),] = as.integer((as.integer(temp1[(temp1$A==(col3-1)&temp1$B==i),]) + as.integer(temp1[(temp1$A==(col3+1)&temp1$B==i),])) / 2)
  temp1[(temp1$A==col4&temp1$B==i),] = as.integer((as.integer(temp1[(temp1$A==(col4-1)&temp1$B==i),]) + as.integer(temp1[(temp1$A==(col4+1)&temp1$B==i),])) / 2)
  
}


temp1 = data.frame(X=paste0(temp1$A, "x", temp1$B), temp1)

temp1$A = NULL
temp1$B = NULL

location <- read.table("/gpfs/ycga/work/liu_yang/dw779/3_FFPE_Project/Image/BLADDER312_position.txt", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])
x = x[-1]

data_filtered <- temp1[temp1$X %in% x,]
write.table(data_filtered, file = 'Filtered_matrix_correct.tsv', sep = '\t',col.names=TRUE, row.names = FALSE,quote = FALSE)

my_data <- read.table(file = 'Filtered_matrix_correct.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
data_filtered <- my_data

#data_filtered <- my_data
##remove Proteins has less than 10 expression
count <- rowSums(data_filtered[,2:ncol(data_filtered)])
data_filtered_binary <- data_filtered[,2:ncol(data_filtered)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)

#log_count <- log(count)

##UMI Count
#region <- 6000  #change the x axis maxium
#test <- data_filtered %>% separate(X, c("A", "B"),  sep = "x")
df <- data.frame(number=1, c=count)
write.csv(df, file="Bladder312_ADT_UMI.csv")


data1 <- read.table("Filtered_matrix_correct.tsv", header = TRUE, sep = "\t", row.names = 1)
data3 <- data.frame(X = rownames(data1), data1)
temp1 <- data3 %>% separate(X, c("A", "B"),  sep = "x")

#repair the strips of the data, here is the x= 34 and 46, y = 20. 

temp1$A = NULL
temp1$B = NULL
temp1$unmapped = NULL
data2 <- t(temp1)
sample1.name <- "Bladder312"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
pbmc          <- CreateSeuratObject(matrix1.data, min.cells = 1, project = sample1.name)

alphabet = c(
  '0' = "#73A1BE",
  '1' = "#EB545C",
  '2' = "#efd510",
  '3' = '#7F3C8D',
  '4' = '#80FF08',
  '5' = '#DBA091')

pbmc <- NormalizeData(pbmc, normalization.method = 'CLR', margin = 2) 
pbmc <- ScaleData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 150)
#pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunPCA(
  pbmc,
  #reduction.name = "apca",
  npcs = 10,
  approx = FALSE,
  verbose = FALSE
)
ElbowPlot(pbmc)

pbmc <- RunUMAP(pbmc, dims = 1:5, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:5, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution=0.2, verbose = FALSE)

p1 <- DimPlot(pbmc, label = TRUE) + NoLegend() + 
  scale_color_manual(values = alphabet[1:(pbmc@active.ident %>% unique %>% length )])
p1
ggsave("ADT312_UMAP_preview.pdf", plot = p1, dpi = 600, bg = NULL, height = 6, width = 8)

ident <- Idents(pbmc)
df <- data.frame(ident[])
df1 <-data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

pdf(file = paste("Preview_Protein-CLR_5_5_0.2—nfeatures150.pdf",sep =""), width=11.6, height=11.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) + 
  scale_color_manual(values = alphabet[1:(pbmc@active.ident %>% unique %>% length )]) + 
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 16, size = 4.5)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  coord_equal(xlim=c(0,73),ylim=c(73,1)) +
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        #axis.text=element_text(size=20),
        #axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),
        #axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +   
  theme(plot.background = element_rect(fill = "black"))
dev.off()

save(pbmc, file="ADT312_pbmc_seuratV5_V1.rdata")

save.image("my_workspace_251227.RData")

###single——protein——plot
#genedata <- pbmc@assays[["RNA"]]@data
genedata <- GetAssayData(pbmc, assay = "RNA", slot = "data")
genedata <- t(genedata)
gene <-as.data.frame(as.matrix(genedata))
gene$X = row.names(gene)
gene$X
test <- gene %>% separate(X, c("A", "B"),  sep = "x")
test$A
write.csv(gene,file="Bladder312_ADT_CLR_expression.csv")


dir.create(paste(data_dir, "/single_protein_CLR",sep =""), showWarnings = FALSE)

setwd(paste(data_dir, "/single_protein_CLR",sep =""))
protein = colnames(test)


for (i in 1:(length(colnames(test))-2)) {
  # test1 = data.frame(A = test$A, B = test$B, C =test[,protein[i]])
  test2 = test[test[protein[i]]>0, ]
  
  p <- ggplot(test2, aes(x = as.numeric(test2$A), y = as.numeric(test2$B), colour=test2[, protein[i]])) +
    #scale_color_gradientn(colours = c("black", "green")) + 
    #scale_color_gradientn(colours = c("blue","green", "red"),
    #oob = scales::squish) +
    scale_color_gradientn(colours = jdb_palette("brewer_spectra"),
                          oob = scales::squish) +
    ggtitle(protein[i]) +
    #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
    geom_point(shape = 16, size = 4)+
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,73),ylim=c(73,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text=element_text(colour="black",size=30),
          axis.title=element_text(colour="black",size=30,face="bold"),
          legend.text=element_text(colour="black",size=30),
          legend.title = element_blank(),
          #legend.title = element_text(colour="black", size=15, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(size =1, fill = NA))
  pdf(file = paste(protein[i], ".pdf",sep =""), width=11.6, height=11.6)
  print(p)
  dev.off()
}


setwd(dir)
save.image("my_workspace.RData")
