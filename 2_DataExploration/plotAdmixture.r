library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyr)
library(tidyverse)

setwd("/labgenomaarea2/saguis-lucero/2_DataExploration/")

### species TI
samplelist <- read_tsv("/labgenomaarea2/saguis-lucero/2_DataExploration/SI_popfile.txt",col_names = "sample")
all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:6){
  data <- read_delim(paste0("TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.",k,".Q"),
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
  scale_fill_manual(values = c("#A94442"))
ggsave("K1_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E"))
ggsave("K2_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 3) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B"))
ggsave("K3_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 4) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B","#5B8063"))
ggsave("K4_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 5) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B","#5B8063","#2E3B4E"))
ggsave("K5_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 6) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B","#5B8063","#2E3B4E","#A94442"))
ggsave("K6_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

### species LW
samplelist <- read_tsv("/labgenomaarea2/saguis-lucero/2_DataExploration/LW_popfile.txt",col_names = "sample")
all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:6){
  data <- read_delim(paste0("LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.",k,".Q"),
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
  scale_fill_manual(values = c("#A94442"))
ggsave("K1_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E"))
ggsave("K2_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 3) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B"))
ggsave("K3_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 4) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B","#5B8063"))
ggsave("K4_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 5) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B","#5B8063","#2E3B4E"))
ggsave("K5_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)

all_data %>%
  filter(k == 6) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width = 1) +
  xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#A94442","#8B6F9E","#D5A84B","#5B8063","#2E3B4E","#A94442"))
ggsave("K6_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.png", width = 14, height = 5)
