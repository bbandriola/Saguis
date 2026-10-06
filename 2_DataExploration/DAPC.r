###############################################################################
# DAPC with adegenet
# Usage:  Rscript dapc_168k_snps_v2.R --vcf snps.vcf.gz --pop populations.txt \
#                                     --max_K 8 [--out dapc_results]
###############################################################################
suppressPackageStartupMessages({
  library(adegenet)
  library(vcfR)
  library(ggplot2)
  library(parallel)
})
set.seed(42)

# Settings
# Command-line arguments (base R, no extra packages needed)
#   --vcf    input VCF (.vcf or .vcf.gz)                 [required]
#   --pop    population file: 2 columns, no header, "sample  group"  [required]
#   --max_K  maximum K to test when estimating K        [default: 8]
#   --out    output directory                           [default: dapc_results]
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
usage <- paste("Usage: Rscript dapc_168k_snps_v2.R --vcf snps.vcf.gz",
               "--pop populations.txt [--max_K 8] [--out dapc_results]")
if ("--help" %in% args) { cat(usage, "\n"); quit(status = 0) }

vcf_file <- get_arg("--vcf")
pop_file <- get_arg("--pop")
max_K    <- suppressWarnings(as.integer(get_arg("--max_K", 8)))
out_dir  <- get_arg("--out", "dapc_results")

if (is.null(vcf_file) || is.null(pop_file)) stop("Missing arguments.\n", usage)
if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)
if (!file.exists(pop_file)) stop("Population file not found: ", pop_file)
if (is.na(max_K) || max_K < 2) stop("--max_K must be an integer >= 2")
ploidy_n     <- 2

max_missing  <- 0.20    # drop SNPs with >20% missing genotypes
min_maf      <- 0.05    # drop SNPs with MAF < 0.05
#max_ind_miss <- 0.30    # drop individuals with >30% missing

estimate_K   <- TRUE    # also estimate K automatically (as a check vs prior groups)
# max_K is set from --max_K above (keep small when you have few individuals)

# Use the cores the scheduler gave you (SLURM), else all-but-one
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
vcf <- read.vcfR(vcf_file, verbose = FALSE)
gl  <- vcfR2genlight(vcf)
ploidy(gl) <- ploidy_n
rm(vcf); invisible(gc())
cat("Loaded:", nInd(gl), "individuals,", nLoc(gl), "SNPs\n")

# Set prior groups 
pops <- read.table(pop_file, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("ind", "pop"))

missing_in_vcf <- setdiff(pops$ind, indNames(gl))
missing_in_pop <- setdiff(indNames(gl), pops$ind)
if (length(missing_in_vcf) || length(missing_in_pop)) {
  stop("Sample name mismatch.\n  In pop file but not VCF: ",
       paste(missing_in_vcf, collapse = ", "),
       "\n  In VCF but not pop file: ",
       paste(missing_in_pop, collapse = ", "))
}
pops <- pops[match(indNames(gl), pops$ind), ]
pop(gl) <- factor(pops$pop)

# Filtering (if necessary)
x0       <- as.matrix(gl)
ind_miss <- rowMeans(is.na(x0)); rm(x0)
#gl       <- gl[ind_miss <= max_ind_miss, ]
gl       <- gl[]
cat("After individual filter:", nInd(gl), "individuals\n")

snp_miss <- glNA(gl) / nInd(gl)
freq     <- glMean(gl)
maf      <- pmin(freq, 1 - freq)
gl       <- gl[, which(snp_miss <= max_missing & maf >= min_maf)]
cat("After SNP filter:", nLoc(gl), "SNPs\n")

grp <- pop(gl)
K   <- nlevels(grp)
N   <- nInd(gl)
cat("Prior groups (K =", K, "):\n"); print(table(grp))

# PCA overview
pca <- glPca(gl, nf = min(N - 1, 10), parallel = TRUE, n.cores = n_cores)
eig_pct <- 100 * pca$eig / sum(pca$eig)

pc_df <- data.frame(ind = indNames(gl), PC1 = pca$scores[, 1],
                    PC2 = pca$scores[, 2], group = grp)
p_pca <- ggplot(pc_df, aes(PC1, PC2, colour = group)) +
  geom_point(size = 3) +
  labs(x = sprintf("PC1 (%.1f%%)", eig_pct[1]),
       y = sprintf("PC2 (%.1f%%)", eig_pct[2]), title = "PCA") +
  theme_bw()
ggsave(file.path(out_dir, "pca.pdf"), p_pca, width = 7, height = 5)

# (Optional) Automatic estimation of K
# choose.n.clust = FALSE + explicit n.pca => runs without any interaction.
if (estimate_K) {
  n_pca_fc <- N - 1
  fc <- find.clusters(gl, max.n.clust = min(max_K, N - 1), n.pca = n_pca_fc,
                      choose.n.clust = FALSE, criterion = "diffNgroup",
                      n.start = 100)
  K_est <- length(unique(fc$grp))
  cat("\nAutomatically estimated K =", K_est, " (prior K =", K, ")\n")
  cat("Inferred clusters vs prior groups:\n")
  print(table(inferred = fc$grp, prior = grp))

  bic <- data.frame(K = seq_along(fc$Kstat), BIC = as.numeric(fc$Kstat))
  write.csv(bic, file.path(out_dir, "find_clusters_BIC.csv"), row.names = FALSE)
  pdf(file.path(out_dir, "find_clusters_BIC.pdf"), width = 6, height = 4.5)
  plot(bic$K, bic$BIC, type = "b", pch = 19, xlab = "Number of clusters (K)",
       ylab = "BIC", main = "BIC vs K")
  abline(v = K_est, lty = 2, col = "red")
  dev.off()
}

# Choose number of PCs by cross-validation
# With N = 15 (5 per group) CV is fragile: training.set = 0.8 keeps ~4 per group
# for training and ~1 per group for validation (stratified by group).
x <- as.matrix(gl)
miss <- which(is.na(x), arr.ind = TRUE)
if (nrow(miss)) {                       # mean-impute missing genotypes
  cm <- colMeans(x, na.rm = TRUE)
  x[miss] <- cm[miss[, 2]]
}

n_pca_max <- min(10, floor(0.8 * N) - 1)

pdf(file.path(out_dir, "xvalDapc.pdf"), width = 7, height = 5)
xv <- xvalDapc(x, grp,
               n.pca.max    = n_pca_max,
               training.set = 0.8,
               result       = "groupMean",
               center       = TRUE, scale = FALSE,
               n.rep        = 50,
               xval.plot    = TRUE,
               parallel     = "multicore", ncpus = n_cores)
dev.off()
rm(x); invisible(gc())

best_npca <- as.integer(xv[["Number of PCs Achieving Lowest MSE"]])[1]
best_npca <- max(best_npca, K - 1)      # need at least K-1 PCs for K-1 DAs
cat("\nNumber of PCs retained:", best_npca, "\n")
print(xv[2:6])

# Final DAPC
dapc_res <- dapc(gl, grp, n.pca = best_npca, n.da = K - 1)

# a-score (checks overfitting; compare optimum to best_npca)
pdf(file.path(out_dir, "a_score.pdf"), width = 7, height = 5)
ascore <- optim.a.score(dapc_res, smart = FALSE, n.sim = 50)
dev.off()

# Plots
# maximum 6, add more if needed
cols <- c("#A94442","#8B6F9E","#D5A84B","#5B8063","#3B5F73","#C97A3E")
  

pdf(file.path(out_dir, "dapc_scatter.pdf"), width = 7, height = 6)
scatter(dapc_res, col = cols, pch = 19, cex = 2, scree.da = TRUE,
        scree.pca = TRUE, posi.pca = "bottomleft", posi.da = "topright",
        legend = TRUE, clabel = 0, posi.leg = "topleft")
dev.off()

pdf(file.path(out_dir, "dapc_compoplot.pdf"), width = 9, height = 5)
compoplot(dapc_res, col = cols, show.lab = TRUE, legend = TRUE,
          posi = "bottomleft", cex.names = 0.7)
dev.off()

# Outputs
write.csv(data.frame(ind = indNames(gl), group = grp, dapc_res$ind.coord),
          file.path(out_dir, "dapc_ind_coords.csv"), row.names = FALSE)

write.csv(data.frame(ind = indNames(gl), original = grp,
                     assigned = dapc_res$assign, round(dapc_res$posterior, 4)),
          file.path(out_dir, "dapc_membership_probs.csv"), row.names = FALSE)

cat("\nOriginal vs assigned group:\n")
print(table(original = grp, assigned = dapc_res$assign))

contrib <- data.frame(snp = locNames(gl), dapc_res$var.contr)
write.csv(contrib, file.path(out_dir, "dapc_snp_contributions.csv"), row.names = FALSE)

for (k in seq_len(K - 1)) {
  thr <- quantile(dapc_res$var.contr[, k], 0.999)
  writeLines(locNames(gl)[dapc_res$var.contr[, k] >= thr],
             file.path(out_dir, sprintf("top_snps_DA%d.txt", k)))
  pdf(file.path(out_dir, sprintf("dapc_loadings_DA%d.pdf", k)), width = 10, height = 4)
  loadingplot(dapc_res$var.contr[, k], threshold = thr,
              main = sprintf("SNP loadings - DA%d", k),
              xlab = "SNP index", ylab = "Contribution", lab.jitter = 0)
  dev.off()
}

saveRDS(dapc_res, file.path(out_dir, "dapc_object.rds"))
sessionInfo()
