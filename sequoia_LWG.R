# Pedigree reconstruction
# Example for LWG samples

# Install packages if not already installed 
if (!requireNamespace("sequoia", quietly = TRUE)) install.packages("sequoia", repos = "https://cloud.r-project.org")
if (!requireNamespace("vcfR",    quietly = TRUE)) install.packages("vcfR",    repos = "https://cloud.r-project.org")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr",   repos = "https://cloud.r-project.org")
if (!requireNamespace("readr",   quietly = TRUE)) install.packages("readr",   repos = "https://cloud.r-project.org")

library(sequoia)
library(vcfR)
library(dplyr)
library(readr)

# Read VCF
vcf <- read.vcfR("1SNPevery25kb_minDP18_LWGsamples.recode.vcf")

# Extract genotypes and convert to numeric (0/1/2 = ref hom/het/alt hom)
# extract.gt() pulls the GT field as character strings ("0/0", "0/1", "1/1").
# The case_when() logic then recodes these to integers: 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate, -9 = missing.
# Since VCF stores SNPs as rows and samples as columns, t() transposes the result so you get samples × SNPs — the format sequoia expects.
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Recode genotypes to counts of alt allele
recode_gt <- function(x) {
  case_when(
    x %in% c("0/0", "0|0") ~ 0L,
    x %in% c("0/1", "1/0", "0|1", "1|0") ~ 1L,
    x %in% c("1/1", "1|1") ~ 2L,
    TRUE ~ -9L   # missing
  )
}

# apply recoding — gt is SNPs x samples, need to transpose to samples x SNPs
gt_num <- apply(gt, 2, function(col) {   # col = one sample
  dplyr::case_when(
    col %in% c("0/0", "0|0") ~ 0L,
    col %in% c("0/1", "1/0", "0|1", "1|0") ~ 1L,
    col %in% c("1/1", "1|1") ~ 2L,
    TRUE ~ -9L
  )
})

# result: SNPs x samples — transpose so rows = individuals
GenoM <- t(gt_num)   # now: samples x SNPs

# Clean up SNP column names (CHROM_POS)
colnames(GenoM) <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")

# Check sample IDs look right
head(rownames(GenoM))

# Optional: export genotype table
sequoia_df <- data.frame(ID = rownames(GenoM), GenoM, check.names = FALSE)
write.csv(sequoia_df,
          file    = "sequoia_vcf_13Ksnps_LWG.csv",
          row.names = FALSE,
          quote   = FALSE)

# Genotype QC
# CheckGeno() examines the matrix and reports: proportion of missing data per individual
# and per SNP, minor allele frequencies, SNPs with more than 2 alleles, and heterozygosity.
# Plot = TRUE produces diagnostic histograms.
CheckGeno(GenoM, Plot = TRUE)

# Load life history table
# Needs columns: ID, Sex (1=F,2=M,3=unk), BirthYear (or NA)
life <- read_delim("lifehistory_LWG.csv", delim = "\t")
life <- as.data.frame(life)
life$BirthYear
which(!life$BirthYear %in% as.integer(life$BirthYear) & !is.na(life$BirthYear))

# Parent assignment only
# Runs sequoia in parentage-assignment mode only: for each individual, it searches
# among all sampled individuals for a parent-offspring pair with a sufficiently high likelihood ratio.
# This is faster than the full run and is a good first check. Err = 0.05 sets a relatively high assumed
# genotyping error rate — reasonable for WGS or RADseq data; lower it if your SNP calls are high-confidence.
ParOUT <- sequoia(
  GenoM        = GenoM,
  LifeHistData = life,
  Module       = "par",
  Err          = 0.05
)

# Full pedigree reconstruction
# Extends parentage assignment with sibship clustering and grandparent assignment.
# For pairs where the parent is not genotyped, sequoia creates "dummy" parent nodes (e.g. F0001, M0001)
# and links sibships through them. This produces a multigenerational pedigree even when not all founders are sampled
SeqOUT <- sequoia(
  GenoM        = GenoM,
  LifeHistData = life,
  Module       = "ped",
  Err          = 0.05
)

# Outputs
# SeqOUT$Pedigree — full reconstructed pedigree (ID, sire, dam for each individual)
# ParOUT$PedigreePar — parentage-only results for comparison
# SeqOUT$MaybeParent — pairs that scored above the filter threshold but below the assignment threshold (ambiguous cases worth inspecting)
# GetMaybeRel() — recalculates likelihood ratios for those ambiguous pairs, helping you decide if they're parent-offspring, siblings, or something else
ped     <- SeqOUT$Pedigree
ped_par <- ParOUT$PedigreePar
may     <- SeqOUT$MaybeParent

GetMaybeRel(
  GenoM        = GenoM,
  SeqList      = SeqOUT,
  LifeHistData = life,
  Err          = 0.05
)

write.csv(SeqOUT$Pedigree,
          file      = "pedigree_assigned_13Ksnps_LWG.csv",
          row.names = FALSE)

