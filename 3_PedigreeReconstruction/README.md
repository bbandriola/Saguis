# Pedigree reconstruction

This document describes the workflow used to reconstruct the pedigree of populations of _Leontocebus weddelli_ and _Tamarinus imperator_ from genome-wide SNP data, using the R package sequoia (doi:10.1111/1755-0998.12665).

### Input data
- VCF file produced from '1_Mapping_SNPCalling'
- Individual metadata including: ID, sex, birth year estimation 

### 1. Create new VCF files for each species
In this case, we are using VCF with the filters listed below and retaining only individuals of a single species (LWGsamples.txt).
- Biallelic sites only
- Minimum genotype quality Q30
- No missing genotypes
- Minor allele frequency (MAF) ≥ 0.05
- Normalized variant representation
- Minimum-depth 15
- Thin 1 SNP every 50 kb

```
vcftools --gzvcf 1SNPevery50kb_minDP15_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.recode.vcf.gz --keep LWGsamples.txt --recode --recode-INFO-all --out 1SNPevery50kb_minDP15_LWGsamples 
```

### 2. Create the metadata file
ID, Sex (1=F,2=M,3=unk), BirthYear (or NA), and optionally BY.min/BY.max (the earliest/latest possible birth year for an individual).

Example:

```
ID      Sex     BirthYear       BY.min  BY.max
LWG1TA1 2       NA      NA      2023
LWG1TA2 2       NA      2026    2026
LWG1TA3 2       NA      2026    2026
LWG1TA4 1       NA      NA      2023
LWG1TA5 1       NA      NA      2023
LWG1TA6 1       NA      2023    2024
LWG2TA1 2       NA      2024    2025
LWG2TA2 1       NA      NA      2023
LWG2TA3 2       NA      NA      2023
LWG3TA1 2       NA      NA      2023
LWG3TA2 1       NA      NA      2023
LWG3TA3 2       NA      2025    2026
LWG3TA4 3       NA      2024    2025
LWG3TA5 2       NA      NA      2023
LWG3TA6 2       NA      2023    2024
LWG3TA7 1       NA      2023    2024
LWG3TA8 2       NA      NA      2023
LWG4TA1 1       NA      2023    2024
LWG4TA2 1       NA      NA      2023
LWG4TA3 2       NA      NA      2023
LWG4TA4 2       NA      2023    2024
LWG4TA5 2       NA      2024    2025
```


### 3. Running sequoia

The full script: _sequoia_LWG.R_

If you run on a cluster
```
Rscript sequoia_LWG.R
```
