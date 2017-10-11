library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)


### Barplot for TF count in each ideas bin:
# weights <- ifelse(pcaOutput2$PC1 < -5 & abs(pcaOutput2$PC2) > 10, 2, 1)

input_file <- "~/Dropbox/encode_3/tf_factor_count_with_ideas_bins/files/final_tf_ideas_piechart_combined_barplot_data.bed"
tf_count_ideas_df <- fread(input_file, sep="\t", header= TRUE)

tf_count_ideas_df <- tf_count_ideas_df %>% 
						data.frame %>%
						arrange(tf_counts)

output_file_name <- paste0("~/Dropbox", "/", "final_tf_count_with_ideas_bin_barplot_1.pdf")				
pdf(output_file_name)
# test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))
reordered_ideas_state = factor(tf_count_ideas_df$ideas_state, levels=tf_count_ideas_df$ideas_state)

barplot <- ggplot(tf_count_ideas_df, aes(x=reordered_ideas_state, y=tf_counts)) + 
	geom_bar(stat="identity", fill = "grey") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    ylab("Transcription Factor Counts") + xlab("IDEAS states") + 
    theme_bw() + 
    ggtitle("Total TF counts across each IDEAS states") + 
    theme(
    axis.text.y = element_text(size=6, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +
    #geom_text(aes(label=tf_counts), color="black", size=2.5, hjust = -0.1 ) +
	theme(axis.text.x = element_text(size=8, angle = 90, vjust = 0.5, hjust = 1)) +
	coord_flip() + 
	#scale_fill_manual(values=custom_col )+
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() 

print(barplot)
dev.off()



# Heatmap:
input_file <- "~/Dropbox/encode_3/tf_factor_count_with_ideas_bins/files/final_tf_ideas_piechart_combined_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
#read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "final_tf_count_with_ideas_bin_heatmap.pdf")				
pdf(output_file_name)

grey_col <-  colorRampPalette("grey")
red_col <- colorRampPalette(c("indianred1"))

custom_col <- c(grey_col(1), red_col(1))

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.7)
heatmap.2(mat_data,
  #dendrogram = "none",
  # Rowv = TRUE, 
  # Colv = TRUE,
  #scale="none",
  main = "TF occupancy across different IDEAS states", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  na.color="grey"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    

dev.off()


