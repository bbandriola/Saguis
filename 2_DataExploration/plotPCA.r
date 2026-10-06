library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyr)

setwd("./")
base_name <- "SI_1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac"

pca <- read_table(paste0(base_name, ".eigenvec"), col_names = FALSE)
eigenval <- scan(paste0(base_name, ".eigenval"))

names(pca)[1:ncol(pca)] <- paste0("Ind")
names(pca)[2:ncol(pca)] <- paste0("Lineages")
names(pca)[1:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# # names(pca)[2] <- paste0("Lineages")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# percentage of variance
pdf(paste0(base_name, ".PCvariance.pdf"),width = 10,height = 6)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
dev.off()
# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

###############################################
########## with individual identification #####
###############################################
pdf(paste0(base_name,".PCA_indID.pdf"),width = 10,height = 10)
b <- ggplot(pca, aes(PC1, PC2, color = Lineages, label = Ind)) + 
  geom_text(position = position_jitter(seed = 1, width = 0.02, height = 0.03), 
            size = 4, show.legend = TRUE) +
  scale_color_manual(values = c("#A94442","#8B6F9E","#D5A84B","#5B8063")) +
  coord_equal() + 
  theme_bw(base_size = 20) +
  theme(axis.ticks = element_line(colour = "white", size = 0.3),
        panel.grid.major = element_line(colour = "white", size = 0.3),
        panel.grid.minor = element_blank())

b + xlab(paste0("PC1(", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2(", signif(pve$pve[2], 3), "%)")) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20))
dev.off()
