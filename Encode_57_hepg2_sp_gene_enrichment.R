library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)
library(xlsx)
library("cluster")

args <-  commandArgs(TRUE)
input_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]

read_file <- fread("~/Dropbox/encode_3/hepg2_liver_gene_sp_analysis/files/final_tf_gene_enrichment_pval_corr_heatmap_data.bed", sep="\t", header=TRUE)
read_df <- as.data.frame(read_file)

#read_df <- read_df[,select_cols]
read_df <- as.data.frame(read_df)


#read_df[is.na(read_df)] <- -0.05

data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data) 


output_file_name <- paste0("~/Dropbox", "/", "hepg2_liver_57_sp_gene_heatmap.pdf")        
pdf(output_file_name)

summary(read_df) #max_val = 325
#pair_break =  c(c(-1,0.099), seq(0,1,length=9))
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.50)
heatmap.2(mat_data,
  main = "TF enrichment over highly liver specific & hepg2 expressed genes (-log2(pvalue))", 
  xlab = "Liver-HepG2 specific genes",
  ylab = "Transcription Factors",
  col=greenred(10),  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #symbreaks = min(mat_data, na.rm=TRUE),
  na.color="grey"
#  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/9))
  )    

dev.off()


###################
# For 30 filtered genes:

read_file <- fread("~/Dropbox/encode_3/hepg2_liver_gene_sp_analysis/files/final_tf_enrichment_pval_corr_heatmap_data_gene_filtered.bed", sep="\t", header=TRUE)
read_df <- as.data.frame(read_file)

#read_df <- read_df[,select_cols]
read_df <- as.data.frame(read_df)


#read_df[is.na(read_df)] <- -0.05

data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data) 


output_file_name <- paste0("~/Dropbox", "/", "hepg2_liver_30_sp_filtered_gene_heatmap.pdf")        
pdf(output_file_name)

summary(read_df) #max_val = 325
#pair_break =  c(c(-1,0.099), seq(0,1,length=9))
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.50)
heatmap.2(mat_data,
  main = "TF enrichment over highly liver specific & hepg2 expressed genes(-log2(pvalue))", 
  xlab = "Liver-HepG2 specific genes",
  ylab = "Transcription Factors",
  col=greenred(10),  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #symbreaks = min(mat_data, na.rm=TRUE),
  na.color="grey"
#  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/9))
  )    

dev.off()


