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
#library(dittoSeq)
library(GenomicRanges)
library(BiocParallel)
library(aCGH)

setwd("")

load("/gpfs/gibbs/project/liu_yang/xd97/CNV_ref_genome/copykit_environment_human.RData")

#prepare RNA data
load("RNAobj_proteinCoding_addCNVclone_V2.rdata")
pbmc=obj_rna


#prepare DNA data
load("./copykit_addRNASeuratClusters_V3.rdata")

###==========================================================================###
###step4.1 gene dosage-calculate mean of RNA expression
###==========================================================================###
head(pbmc@meta.data)
#save(pbmc, file = "RNA552_pbmc_celltype_addCNVsubclone_20260117.rdata")

#remove NA in subclones
#pbmc <- pbmc[, !is.na(pbmc@meta.data$subclones)]
#save(pbmc, file = "RNA552_pbmc_celltype_addrevisedCNVsubclone_rmNA_20260117.rdata")

#revised
head(pbmc@meta.data)

#rna_meta_revised <- data.frame(pbmc@meta.data)
#head(rna_meta_revised)

#rna_meta_revised <- rna_meta_revised %>%
#  dplyr::mutate(
#    subclones = ifelse(
#      celltype %in% c("Epithelial cells", "Muscle cells"),
#      "normal cells",
#      "c1"
#    )
#  )

#all(rownames(rna_meta_revised) == colnames(pbmc))
#rna_meta_revised <- rna_meta_revised[colnames(pbmc), ]
#pbmc@meta.data <- rna_meta_revised
#head(pbmc@meta.data)
###===============================================================###
#log normalize RNA data for compare
pbmc <- NormalizeData(pbmc, assay = "RNA")

subclones <- as.vector(pbmc@meta.data$subclones)
# Gene expression matrix for each cell,必须是RNA assay
object_bulk <- as.data.frame(t(GetAssayData(pbmc, assay = "RNA", slot = "data")))

# Add subclone index
object_bulk$subclone <- subclones

#calculate the mean expression of RNA in eaxh subclone
# Replace zeros with NA for mean and median calculations
all_genes <- colnames(object_bulk)
# mean expression of each subclone
# object_bulk_remove_zero <- object_bulk %>%
#       dplyr::mutate(across(where(is.numeric), ~replace(., . == 0, NA)))
object_bulk_mean <- object_bulk %>%
  dplyr:: group_by(subclone) %>%
  dplyr:: summarise(across(everything(),   ~mean(., na.rm = TRUE)))
object_bulk_mean <- as.data.frame(t(object_bulk_mean))
# Set column names and remove first row (which contains subclone names)
#object_bulk_mean <- object_bulk_mean %>% column_to_rownames("subclones")
colnames(object_bulk_mean) <- object_bulk_mean[1, ]
object_bulk_mean <- object_bulk_mean[-1, ]
object_bulk_mean$gene <- rownames(object_bulk_mean)
head(object_bulk_mean)
# Get rid of unassigned subclone cells
#object_bulk_mean <- object_bulk_mean %>%
#  select(-matches("^(unassigned|unsigned)$"))

###==========================================================================###
###step4.2 gene dosage-process DNA data
###==========================================================================###
copykit
varbin_mtx_all_log2 <- copykit

varbin_mtx_all_log2  <- calcConsensus(varbin_mtx_all_log2 )
varbin_mtx_all_log2  <- runConsensusPhylo(varbin_mtx_all_log2 )

raw_mtx_cnv <- as.data.frame(varbin_mtx_all_log2@consensus)
#add bin name
chr_ranges <-as.data.frame(SummarizedExperiment::rowRanges(varbin_mtx_all_log2))
raw_mtx_cnv$bin <- chr_ranges$seqnames
raw_mtx_cnv$abspos <- chr_ranges$abspos
head(raw_mtx_cnv)


raw_mtx_cnv <- raw_mtx_cnv %>%
  arrange(bin, abspos) %>%                     
  mutate(bin = paste0(bin, ".", row_number()))  

raw_mtx_cnv <- as.data.frame(raw_mtx_cnv)
rownames(raw_mtx_cnv) <- raw_mtx_cnv$bin
raw_mtx_cnv$bin <- NULL
head(raw_mtx_cnv)

#colnames(raw_mtx_cnv)[colnames(raw_mtx_cnv) == "normal cells"] <- "Diploid"

#write.csv(raw_mtx_cnv, 
#          file = "mouseliver552_raw_mtx_cnv_concensus-260117.csv", 
#          row.names = FALSE)


cnv_changes_diploid <- raw_mtx_cnv

clone_cols <- setdiff(colnames(raw_mtx_cnv), c("C1", "abspos","Diploid"))
clone_cols

for (clone in clone_cols) {
  new_col_name <- paste0("cnv_", clone, "_vsDiploid")
  cnv_changes_diploid[[new_col_name]] <- raw_mtx_cnv[[clone]] - raw_mtx_cnv$C1
}


#Plot diff_CNV
cnv_changes_diploid

###==========================================================================###
###step4.3 build gene-bin map, only run once for each ref-mm10,hg19,or hg38
###==========================================================================###
gene_bin <- read.table("/gpfs/ycga/work/liu_yang/xd97/FFPE_triomics_analysis/1_FFPE_CNV_final/14_BRX50_PBRD_BM01_BM02_integration/step3_mRNA_genedosage_7clones/build_hg38_gene_bin_map/hg38_bin_gene_map_geneStart.tsv", header = 1)


###==========================================================================###
###step4.4 add coorelation between CMV bins and gene based on bin_gene_map3
###==========================================================================###
head(cnv_changes_diploid)

#add bins information
cnv_changes_diploid$bins <- as.numeric(gsub(".*\\.", "", rownames(cnv_changes_diploid)))

# find genes for each bin
#gene_bin <- read.table("/gpfs/ycga/work/liu_yang/xd97/FFPE_triomics_analysis/1_FFPE_CNV_final/3_MouseLungT_merged3batch/visualization2_filterbam1st/gene_dosage_mm10_2_plots/build_bin_miRNA_map/mm10_binNumber_protein_coding_map.tsv", header = 1)

gene_bin_sel <- gene_bin %>% 
  dplyr::filter(gene %in% all_genes ) #all_genes from RNA expression data:object_bulk

names(gene_bin_sel) <- c("gene", "bins")
#connect genes and cnv based on bins number information
cnv_genes_mtx_diploid <- left_join(cnv_changes_diploid, gene_bin_sel, by = "bins")

cnv_genes_mtx_diploid
# find expression data for each gene
cnv_genes_mean_exp_mtx_diploid <- left_join(cnv_genes_mtx_diploid,object_bulk_mean, by = "gene" )
names(cnv_genes_mean_exp_mtx_diploid) <- gsub(".y$", "_exp", names(cnv_genes_mean_exp_mtx_diploid))
names(cnv_genes_mean_exp_mtx_diploid) <- gsub(".x$", "_cnv", names(cnv_genes_mean_exp_mtx_diploid))
head(cnv_genes_mean_exp_mtx_diploid)

#colnames(cnv_genes_mean_exp_mtx_diploid)[colnames(cnv_genes_mean_exp_mtx_diploid) == "normal cells"] <- "Diploid_exp"

#write.csv(cnv_genes_mean_exp_mtx_diploid, 
#          file = "mouselung552_cnv_genes_mean_exp_mtx_diploid-260117.csv", 
#          row.names = TRUE)

#MEAN
# calculate mean expression of each bin
exp_columns <- grepl("exp$",colnames(cnv_genes_mean_exp_mtx_diploid))
cnv_genes_mean_exp_mtx_diploid[exp_columns] <- lapply(cnv_genes_mean_exp_mtx_diploid[exp_columns], as.numeric)

cnv_genes_mean_exp_mtx_diploid <- cnv_genes_mean_exp_mtx_diploid %>%
  dplyr::mutate(across(ends_with("_exp"), ~ ifelse(. > 0.05, ., NA))) # set a cutoff to filter low expression genes

gene_per_bin <- cnv_genes_mean_exp_mtx_diploid %>%
  dplyr::group_by(bins) %>%
  dplyr::summarise(count=n()) 

median(gene_per_bin$count)
# calculate expression minus c1 for each clone
clone_exp_cols <- grep("_exp$", colnames(cnv_genes_mean_exp_mtx_diploid), value = TRUE)
clone_exp_cols <- setdiff(clone_exp_cols, c("Diploid_exp", "C1_exp"))
clone_exp_cols

for (col in clone_exp_cols) {
  new_col_name <- paste0("sub_", col) 
  cnv_genes_mean_exp_mtx_diploid[[new_col_name]] <- cnv_genes_mean_exp_mtx_diploid[[col]] - cnv_genes_mean_exp_mtx_diploid$C1_exp
}

head(cnv_genes_mean_exp_mtx_diploid)
#write.csv(cnv_genes_mean_exp_mtx_diploid, 
#          file = "mouselung_cnv_genes_mean_exp_mtx_diploid_RNAdiff—cuoff0.05-260117.csv", 
#          row.names = FALSE)


clone_name <- c("C2","C3","C4","C5","C6","C7")

for (i in 1:length(clone_name)) {
  
  cnv_col <- paste0(clone_name[i], "_cnv")        
  group_col <- paste0("group_", clone_name[i])   
  
  group <- 1 
  cnv_genes_mean_exp_mtx_diploid[[group_col]] <- NA  
  

  for (j in 1:(nrow(cnv_genes_mean_exp_mtx_diploid)-1)) {
    
  
    cnv_genes_mean_exp_mtx_diploid[[group_col]][j] <- group
    
   
    if (cnv_genes_mean_exp_mtx_diploid[[cnv_col]][j+1] == cnv_genes_mean_exp_mtx_diploid[[cnv_col]][j]) {
      cnv_genes_mean_exp_mtx_diploid[[group_col]][j+1] <- group
    } else {
      group <- group + 1   # CNV 变化，group +1
      cnv_genes_mean_exp_mtx_diploid[[group_col]][j+1] <- group
    }
  }
}

head(cnv_genes_mean_exp_mtx_diploid)


sample_mean_mtx_name <- vector("character", length = 0)
clone_name <- c("C2","C3","C4","C5","C6","C7")  
for (i in 1:length(clone_name)) {
  group_col <- paste0("group_", clone_name[i])    
  exp_col <- paste0("sub_", clone_name[i], "_exp") 
  cnv_col <- paste0(clone_name[i], "_cnv")    
  diff_cnv_col <- paste0("cnv_",clone_name[i],"_vsDiploid" ) 
  name <- paste0("cnv_genes_exp_mean_mtx_bin_", clone_name[i])
  
  mtx <- data.frame(
    group = cnv_genes_mean_exp_mtx_diploid[[group_col]],
    cnv = cnv_genes_mean_exp_mtx_diploid[[cnv_col]], #cnv 
    sub_bin = cnv_genes_mean_exp_mtx_diploid[[exp_col]], #RNA diff between clone and diploid
    bins = cnv_genes_mean_exp_mtx_diploid$bins,
    #bin = cnv_genes_mean_exp_mtx_diploid$bin,
    #gene_start = cnv_genes_mean_exp_mtx_diploid$gene_start,
    gene = cnv_genes_mean_exp_mtx_diploid$gene,
    diff_cnv = cnv_genes_mean_exp_mtx_diploid[[diff_cnv_col]]
  )
  
 
  mtx_group <- mtx %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(sub_group = median(sub_bin, na.rm = TRUE))
  
 
  mtx <- left_join(mtx, mtx_group, by = "group")
  
 
  write.csv(mtx, file = paste0(name, ".csv"), row.names = FALSE)
  
  assign(name, mtx)
  sample_mean_mtx_name <- append(sample_mean_mtx_name, name)
}

