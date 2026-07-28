library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyr)
library(tidyverse)

k = c(1:6)
cv1.err = c(1.19489,0.10536,0.11482,0.10965,0.11605,0.12485)

pdf("CVK1-6.pdf",)
plot(k,cv1.err, type="b",xlab = "K", ylab = "cross-validation error", main = "CV test for K=1-20",
     cex.main=1.5, cex.lab=1.5,yaxt="n",xaxt="n")
axis(1,cex.axis=1.3)
axis(2,cex.axis=1.3)
dev.off()

samplelist <- read_tsv("sampleorder.txt",
                       col_names = "sample")

all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:6){
  data <- read_delim(paste0("1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.",k,".Q"),
                  col_names = paste0("Q",seq(1:k)),
                  delim=" ")
  data$sample <- samplelist$sample
  data$k <- k
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}

all_data %>%
  filter(k == 1) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))
ggsave("K1_1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))
ggsave("K2_1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 3) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))
ggsave("K3_1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 4) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))
ggsave("K4_1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 5) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))
ggsave("K5_1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 6) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))
ggsave("K6_1SNPevery10kb_minDP10_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)
