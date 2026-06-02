###Aim: this scrpt aims to generating RNA type-specific expression matrices.
###Author: Xiangjun Di

```r
expmat <- read.table("expmat.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

expmat$Gene <- sub("__.*$", "", expmat[,1])

anyNA(expmat$Gene)

expmat <- expmat[,-1]

library(dplyr)

final_expmat <- expmat %>%
  group_by(Gene) %>%
  summarise(across(everything(), sum))

write.table(final_expmat, "finalexpmat_merged.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

library(dplyr)

expmat <- read.table("finalexpmat_merged.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

genetype_dict <- read.table("path/name2type_hsa_dictionary.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("Gene", "Genetype"))

unique(genetype_dict$Genetype)

expmat <- expmat %>%
  left_join(genetype_dict, by = c("Gene" = "Gene"))

keep_types <- c("lncRNA", "miRNA", "misc_RNA", "protein_coding",
                "rRNA", "scaRNA", "snoRNA", "snRNA", "tRNA",
                "vault_RNA", "Y_RNA", "scRNA", "exon", "pseudogene")

expmat$Genetype <- ifelse(expmat$Genetype %in% keep_types, expmat$Genetype, "other RNA")

for(gt in unique(expmat$Genetype)){
  df_sub <- expmat %>% filter(Genetype == gt) %>% select(-Genetype)
  outfile <- paste0(gt, ".tsv")
  write.table(df_sub, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
}

```
