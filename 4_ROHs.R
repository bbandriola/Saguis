# generated with bcftools roh --rec-rate 1e-8 /labgenomaarea2/saguis-lucero/1_Mapping_SNPCalling/bcftools/NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.vcf.gz -o ROHsperind.txt
# Load libraries and read data
library(tidyverse)
#install.packages("ggrepel")
library(ggrepel)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 2_ROHsPlot.r <FilteredPCAandUCE_Max30missinessDepthmin10max100_Sspeluncaecomplexedited.roh.population.txt>")
}

input_file <- args[1]
# Read data with read_delim() for better control over input file parsing
roh <- read_delim(input_file, delim = "\t", col_names = c("Sample", "Chromosome", "RoH_length", "Population"))

# Compute NROH and SROH
roh <- roh %>%
  mutate(RoH_length = as.numeric(RoH_length))

nroh <- roh %>%
  group_by(Sample) %>%
  summarize(NROH = n())

sroh <- roh %>%
  group_by(Sample) %>%
  summarize(SROH = sum(RoH_length))

# Combine population information
data <- inner_join(nroh, sroh, by = "Sample") %>%
  inner_join(roh %>% select(Sample, Population) %>% distinct(), by = "Sample")

# Inspect the resulting dataset to ensure it is correct
#head(data)

my.colors <- c("#0660FB","#4D9F8C", "#7570B3", "#E8488B", "#66A61D","#E93F33","#E6AC2F")

# Create the plot using the preprocessed dataset
snroh_plot <- data %>%
  ggplot(aes(x = SROH, y = NROH, color = Population)) +
  geom_point(size = 3) +
  scale_colour_manual(values =my.colors)+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "NROH vs. SROH", color = "Population")

# Print the plot
print(snroh_plot)
ggsave("dsroh_nroh.png", snroh_plot, width = 8, height = 6, dpi = 300)

#summed_roh%%R

# Create RoH length categories
#roh_categories <- roh %>%
#  mutate(RoH_category = case_when(
#      RoH_length >= 1000000 ~ "< 1 Mb",
#      RoH_length >= 1000000 & RoH_length < 3000000 ~ "1 Mb - 3 Mb",
#      RoH_length >= 3000000 & RoH_length < 5000000 ~ "3 Mb - 5 Mb",
#      RoH_length >= 5000000 ~ "> 5 Mb",
#      TRUE ~ NA_character_)) %>% # Assign NA to lengths < 1 Mb
#  mutate(RoH_category = factor(RoH_category, levels = c("> 5 Mb", "3 Mb - 5 Mb", "1 Mb - 3 Mb","< 1 Mb"))) %>%
#  filter(!is.na(RoH_category)) # Remove rows where RoH_category is NA
roh_categories <- roh %>%
  mutate(
    RoH_length = as.numeric(gsub("[^0-9.]", "", RoH_length)),
    RoH_category = case_when(
      RoH_length < 1e5 ~ "< 100 kb",
      RoH_length >= 1e5 & RoH_length < 5e5 ~ "100 kb - 500 kb",
      RoH_length >= 5e5 & RoH_length < 1e6 ~ "500 kb - 1 Mb",
      RoH_length >= 1e6 & RoH_length < 2e6 ~ "1 Mb - 2 Mb",
      RoH_length >= 2e6 ~ "> 2 Mb",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    RoH_category = factor(
      RoH_category,
      levels = c("> 2 Mb","1 Mb - 2 Mb","500 kb - 1 Mb","100 kb - 500 kb","< 100 kb")
    )
  ) %>%
  filter(!is.na(RoH_category))


# Summing RoH_length by Sample, RoH_category, and Population
#summed_roh <- roh_categories %>%
#  group_by(Sample, RoH_category, Population) %>%
#  summarise(total_RoH_length = sum(RoH_length), .groups = "drop")  # Sum the RoH_length for each group
summed_roh <- roh_categories %>%
  group_by(Sample, RoH_category, Population) %>%
  summarise(total_RoH_length = sum(RoH_length), .groups = "drop")
my_sample_order <- c("LWG1TA1","LWG1TA2","LWG1TA3","LWG1TA4","LWG1TA5","LWG1TA6","LWG2TA1","LWG2TA2","LWG2TA3","LWG3TA1","LWG3TA2","LWG3TA3","LWG3TA4","LWG3TA5","LWG3TA6","LWG3TA7","LWG3TA8","LWG4TA1","LWG4TA2","LWG4TA3","LWG4TA4","LWG4TA5","SIG1TA1","SIG1TA2","SIG1TA3","SIG1TA4","SIG1TA5","SIG2TA1","SIG2TA2","SIG2TA3","SIG2TA4","SIG2TA5","SIG3TA1","SIG3TA2","SIG3TA3","SIG3TA4","SIG3TA5")  # your list here
my_sample_order_rev <- rev(my_sample_order)
summed_roh <- summed_roh %>%
  filter(Sample %in% my_sample_order) %>%
  mutate(Sample = factor(Sample, levels = my_sample_order_rev))
roh_colors <- c(
  "< 100 kb"      = "#db8a74",
  "100 kb - 500 kb" = "#a7333f",
  "500 kb - 1 Mb"   = "#74121d",
  "1 Mb - 2 Mb"   = "#c52233",
  "> 2 Mb"        = "#370617"
)
print(summed_roh)

# Create the stacked bar plot
stacked_barplot <- summed_roh %>%
  ggplot(aes(x = Sample, y = total_RoH_length/1000000, fill = RoH_category)) +
  geom_bar(stat = "identity",width = 0.8) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(hjust = 1),
    strip.text = element_text(angle = 90, size = 10),panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = roh_colors) +
  labs(
    x = "Sample",
    y = "Total RoH Length (Mb)",
    fill = "RoH Category"
  )+ coord_flip()
# Print the plot
print(stacked_barplot)
ggsave("smallcategories_rohs_barplot.pdf", stacked_barplot, width = 8, height = 6, dpi = 300)
