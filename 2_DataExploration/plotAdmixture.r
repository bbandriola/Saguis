# load tidyverse package
library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyr)

k = c(6:20)
cv1.err = c(0.42939,0.53406,0.56896,0.59632,0.50977,0.47235,0.42158,
            0.60191,0.59212,0.64720,0.78124,0.62473,0.43106,0.62745,0.58790)

pdf("CVK6-20_ALLspecies.pdf",)
plot(k,cv1.err, type="b",xlab = "K", ylab = "cross-validation error", main = "CV test for K=1-20",
     cex.main=1.5, cex.lab=1.5,yaxt="n",xaxt="n")
axis(1,cex.axis=1.3)
axis(2,cex.axis=1.3)
dev.off()

tbl=read.table("1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.1.Q")
tbl$V1 <- factor(tbl$V1, levels = unique(tbl$V1))
tbl_long <- pivot_longer(tbl, cols = V2, names_to = "K", values_to = "Proportion")
p<-ggplot(tbl_long, aes(x = V1, y = Proportion, fill = K)) +
  geom_col(width = 1, color = NA) +  # Full-width bars
  scale_fill_manual(
    values = c("V2" = "#66A61D"),
    labels = c("K=1")
  ) +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave("K1_SspeluncaeComplex.png", plot = p, width = 14, height = 5)

tbl=read.table("1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.2.Q")
tbl$V1 <- factor(tbl$V1, levels = unique(tbl$V1))
tbl_long <- pivot_longer(tbl, cols = V2:V3, names_to = "K", values_to = "Proportion")
p<-ggplot(tbl_long, aes(x = V1, y = Proportion, fill = K)) +
  geom_col(width = 1, color = NA) +  # Full-width bars
  scale_fill_manual(
    values = c("V2" = "#0660FB","V3" = "#4D9F8C"),
    labels = c("K=1","K=2")
  ) +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave("K2_SspeluncaeComplex.png", plot = p, width = 14, height = 5)


tbl=read.table("1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.3.Q")
tbl$V1 <- factor(tbl$V1, levels = unique(tbl$V1))
tbl_long <- pivot_longer(tbl, cols = V2:V4, names_to = "K", values_to = "Proportion")
p<-ggplot(tbl_long, aes(x = V1, y = Proportion, fill = K)) +
  geom_col(width = 1, color = NA) +  # Full-width bars
  scale_fill_manual(
    values = c("V2" = "#0660FB","V3" = "#4D9F8C","V4" = "#7570B3"),
    labels = c("K=1","K=2","K=3")
  ) +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave("K3_SspeluncaeComplex.png", plot = p, width = 14, height = 5)

tbl=read.table("1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.4.Q")
tbl$V1 <- factor(tbl$V1, levels = unique(tbl$V1))
tbl_long <- pivot_longer(tbl, cols = V2:V5, names_to = "K", values_to = "Proportion")
p<-ggplot(tbl_long, aes(x = V1, y = Proportion, fill = K)) +
  geom_col(width = 1, color = NA) +  # Full-width bars
  scale_fill_manual(
    values = c("V2" = "#4D9F8C","V3" = "#0660FB","V4" = "#7570B3","V5" = "deeppink"),
    labels = c("K=1","K=2","K=3","K=4","K=5")
  ) +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave("K4_SspeluncaeComplex.png", plot = p, width = 14, height = 5)

tbl=read.table("1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.5.Q")
tbl$V1 <- factor(tbl$V1, levels = unique(tbl$V1))
tbl_long <- pivot_longer(tbl, cols = V2:V6, names_to = "K", values_to = "Proportion")
p<-ggplot(tbl_long, aes(x = V1, y = Proportion, fill = K)) +
  geom_col(width = 1, color = NA) +  # Full-width bars
  scale_fill_manual(
    values = c("V2" = "#66A61D","V3" = "#0660FB","V4" = "deeppink","V5" = "#7570B3", "V6" = "#4D9F8C"),
    labels = c("K=1","K=2","K=3","K=4","K=5")
  ) +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave("K5_SspeluncaeComplex.png", plot = p, width = 14, height = 5)


tbl=read.table("1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.6.Q")
tbl$V1 <- factor(tbl$V1, levels = unique(tbl$V1))
tbl_long <- pivot_longer(tbl, cols = V2:V7, names_to = "K", values_to = "Proportion")
p<-ggplot(tbl_long, aes(x = V1, y = Proportion, fill = K)) +
  geom_col(width = 1, color = NA) +  # Full-width bars
  scale_fill_manual(
    values = c("V2" = "#66A61D","V3" = "#0660FB","V4" = "deeppink",
    "V5" = "#7570B3", "V6" = "#4D9F8C","V7" = "red"),
    labels = c("K=1","K=2","K=3","K=4","K=5","K=6")
  ) +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave("K6_SspeluncaeComplex.png", plot = p, width = 14, height = 5)


