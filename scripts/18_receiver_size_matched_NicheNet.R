

#!/usr/bin/env Rscript
# SCRIPT 18 — Receiver-size-matched NicheNet sensitivity analysis
# Exploratory sensitivity analysis; does not overwrite frozen analyses.
# Run from the Paper3_Sensory_Islet project root.
# Method: sample 22 genes from the independent 697-gene receiver using
# expression-decile quotas, and compare each subset to ONE shared expression-
# matched null distribution of 22-gene sets from the GSE221156 background.
# The shared null is valid here because all sampled sets use identical quotas.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(nichenetr))

read_positive_int <- function(key, default) {
  x <- suppressWarnings(as.integer(Sys.getenv(key, as.character(default))))
  if (length(x) != 1L || is.na(x) || x < 1L) stop("Invalid ", key)
  x
}

n_subsets <- read_positive_int("N_SUBSETS", 200L)
n_null <- read_positive_int("N_NULL", 2000L)
n_cores <- read_positive_int("N_CORES", 4L)
chunk_size <- read_positive_int("CHUNK_SIZE", 100L)
run_tag <- Sys.getenv("RUN_TAG", "pilot")
if (!grepl("^[A-Za-z0-9_-]+$", run_tag)) stop("RUN_TAG must contain only letters, digits, _ or -")
if (n_null < 999L) warning("Empirical P-value resolution is coarse with fewer than 999 null sets")
seed <- 180221156L
n_bins <- 10L
n_match <- 22L

root <- getwd()
lt_file <- "data/raw/NicheNet/ligand_target_matrix_nsga2r_final.rds"
de_file <- "results/validation/GSE221156/Table_P3_GSE221156_Beta_T2D_vs_Normal_edgeR.csv"
independent_file <- "data/processed/GSE221156_Independent_NicheNet_16a.rds"
broad_file <- "data/processed/BroadLigandUniverse_NicheNet_16b.rds"
original_null_file <- "results/validation/GSE221156/NicheNet_permutation/final/Table_P3_GSE221156_NicheNet_ExpressionMatched_EmpiricalNull.csv"
discovery_null_file <- "results/validation/GSE221156/NicheNet_permutation/final/Table_P3_GSE81608_NicheNet_ExpressionMatched_EmpiricalNull.csv"
outdir <- file.path("results/validation/GSE221156/NicheNet_receiver_size_matched", run_tag)
required <- c(lt_file, de_file, independent_file, broad_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required input files:\n", paste(missing, collapse = "\n"))

expected_md5 <- "b09606b04b2d4490418d9028c0e58b9f"
observed_md5 <- unname(tools::md5sum(lt_file))
if (!identical(observed_md5, expected_md5)) stop("Frozen NicheNet prior MD5 mismatch")
lt <- readRDS(lt_file)
if (!identical(dim(lt), c(33354L, 1226L))) stop("Unexpected frozen NicheNet prior dimensions")

independent <- readRDS(independent_file)
broad <- readRDS(broad_file)
ligands <- sort(unique(as.character(broad$broad_ligands)))
if (length(ligands) != 42L || !all(ligands %in% colnames(lt))) {
  stop("The frozen broad universe does not contain 42 evaluable ligands")
}
receiver <- unique(as.character(independent$receiver_nichenet))
background <- unique(as.character(independent$background_nichenet))
if (length(receiver) != 697L || length(background) != 13712L ||
    !all(receiver %in% background) || !all(background %in% rownames(lt))) {
  stop("Frozen independent receiver/background mismatch")
}

de <- read.csv(de_file, check.names = FALSE)
if (!all(c("gene", "logCPM", "logFC", "FDR") %in% names(de))) stop("DE table missing required columns")
de$gene <- trimws(as.character(de$gene))
de <- de[!is.na(de$gene) & nzchar(de$gene) & !duplicated(de$gene), , drop = FALSE]
if (!setequal(intersect(de$gene, rownames(lt)), background)) {
  stop("GSE221156 background differs from the frozen independent analysis")
}
expr <- de$logCPM
names(expr) <- de$gene
expr <- expr[background]
if (any(!is.finite(expr))) stop("Nonfinite logCPM values in the background")

# Exact decile assignment from original Script 16c (rank + gene-name tie break).
make_expression_bins <- function(background_genes, expression_values, n_bins = 10L) {
  expr <- expression_values[background_genes]
  if (anyNA(expr)) stop("Missing expression value")
  ord <- order(expr, names(expr))
  n <- length(ord)
  bin_ord <- pmin(n_bins, ceiling(seq_len(n) * n_bins / n))
  bins <- integer(n)
  bins[ord] <- bin_ord
  names(bins) <- names(expr)
  bins
}
bins <- make_expression_bins(background, expr, n_bins)
receiver_bin_n <- tabulate(bins[receiver], nbins = n_bins)

# Allocate 22 genes across expression bins in proportion to the 697-gene
# receiver. Largest-remainder rounding guarantees sum(quota) == 22.
ideal <- n_match * receiver_bin_n / length(receiver)
quota <- as.integer(floor(ideal))
remainder <- n_match - sum(quota)
if (remainder > 0L) {
  ord <- order(-(ideal - quota), seq_along(quota))
  quota[ord[seq_len(remainder)]] <- quota[ord[seq_len(remainder)]] + 1L
}
if (sum(quota) != n_match || any(quota > receiver_bin_n)) stop("Invalid bin quota")

sample_stratified <- function(pool, bins, quota) {
  ans <- unlist(lapply(seq_along(quota), function(b) {
    if (quota[b] == 0L) return(character())
    available <- pool[bins[pool] == b]
    if (length(available) < quota[b]) stop("Insufficient genes in expression bin ", b)
    sample(available, quota[b], replace = FALSE)
  }), use.names = FALSE)
  if (length(ans) != sum(quota) || anyDuplicated(ans)) stop("Invalid sampled receiver")
  ans
}

# Restrict matrix once, as in Script 16c, to keep repeated calls efficient.
lt_sub <- lt[background, ligands, drop = FALSE]
rm(lt)
gc(verbose = FALSE)

run_activity <- function(geneset) {
  x <- nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = lt_sub,
    potential_ligands = ligands
  )
  v <- x$aupr_corrected
  names(v) <- as.character(x$test_ligand)
  v <- as.numeric(v[ligands])
  if (length(v) != length(ligands) || any(!is.finite(v))) stop("Invalid ligand activity")
  v
}

# Generate deterministic gene sets in the parent R process; results do not
# depend on the number of parallel workers or the order of task completion.
set.seed(seed)
subset_sets <- replicate(n_subsets, sample_stratified(receiver, bins, quota), simplify = FALSE)
set.seed(seed + 1L)
null_sets <- replicate(n_null, sample_stratified(background, bins, quota), simplify = FALSE)

# Separate checkpoints for nulls and subsets, scoped to this new analysis.
# Input MD5 and all run parameters are included in the checkpoint identity.
fingerprint <- paste(
  observed_md5,
  unname(tools::md5sum(de_file)),
  unname(tools::md5sum(independent_file)),
  unname(tools::md5sum(broad_file)),
  n_subsets, n_null, n_match, seed,
  paste(quota, collapse = ","),
  paste(ligands, collapse = ","),
  sep = "|"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

run_batch <- function(sets, label) {
  checkpoint <- file.path(outdir, paste0("checkpoint_", label, ".rds"))
  scores <- matrix(numeric(), nrow = 0L, ncol = length(ligands),
                   dimnames = list(NULL, ligands))
  if (file.exists(checkpoint)) {
    old <- readRDS(checkpoint)
    if (!identical(old$fingerprint, fingerprint) || !identical(old$label, label)) {
      stop("Existing ", label, " checkpoint belongs to different settings; use a different RUN_TAG")
    }
    scores <- old$scores
    if (ncol(scores) != length(ligands) || !identical(colnames(scores), ligands)) {
      stop("Malformed checkpoint")
    }
    cat("Resuming ", label, " from ", nrow(scores), " / ", length(sets), "\n", sep = "")
  }
  completed <- nrow(scores)
  if (completed > length(sets)) stop("Checkpoint exceeds requested number of sets")
  if (completed < length(sets)) {
    for (start in seq.int(completed + 1L, length(sets), by = chunk_size)) {
      idx <- seq.int(start, min(start + chunk_size - 1L, length(sets)))
      t0 <- proc.time()[["elapsed"]]
      ans <- parallel::mclapply(idx, function(i) run_activity(sets[[i]]),
                                mc.cores = n_cores, mc.preschedule = TRUE,
                                mc.set.seed = FALSE)
      bad <- vapply(ans, function(z) inherits(z, "try-error") ||
                      !is.numeric(z) || length(z) != length(ligands) ||
                      any(!is.finite(z)), logical(1))
      if (any(bad)) stop("NicheNet worker failure in ", label, " chunk ", start)
      chunk <- do.call(rbind, ans)
      colnames(chunk) <- ligands
      scores <- rbind(scores, chunk)
      saveRDS(list(fingerprint = fingerprint, label = label, scores = scores), checkpoint)
      cat(sprintf("%s: %d/%d | chunk %.1f sec\n", label, nrow(scores),
                  length(sets), proc.time()[["elapsed"]] - t0))
    }
  }
  scores
}

cat("\nSCRIPT 18 — RECEIVER-SIZE-MATCHED NICHE-NET SENSITIVITY\n")
cat("Project:", root, "\n")
cat("NicheNet version:", as.character(packageVersion("nichenetr")), "\n")
cat("Frozen prior MD5:", observed_md5, "\n")
cat("Independent receiver:", length(receiver), "genes\n")
cat("Independent background:", length(background), "genes\n")
cat("Size-matched receiver:", n_match, "genes\n")
cat("Expression-bin quotas:", paste(quota, collapse = ", "), "\n")
cat("Ligands:", length(ligands), "\n")
cat("Subsets:", n_subsets, "| null sets:", n_null, "| cores:", n_cores, "\n")
cat("Output:", outdir, "\n\n")

null_scores <- run_batch(null_sets, "expression_matched_null")
subset_scores <- run_batch(subset_sets, "receiver_subsets")

# Empirical one-sided tests per sampled receiver, exactly matching Script 16c:
# P = (1 + number of null AUPR >= observed AUPR)/(1 + N_NULL).
# BH adjustment is performed across the 42 ligands within EACH sampled set.
pvals <- matrix(NA_real_, n_subsets, length(ligands), dimnames = list(NULL, ligands))
fdr <- pvals
ranks <- pvals
for (i in seq_len(n_subsets)) {
  exceed <- vapply(seq_along(ligands), function(j) {
    sum(null_scores[, j] >= subset_scores[i, j])
  }, integer(1))
  pvals[i, ] <- (1 + exceed) / (1 + n_null)
  fdr[i, ] <- p.adjust(pvals[i, ], method = "BH")
  ranks[i, ] <- rank(-subset_scores[i, ], ties.method = "average")
}

# Rank similarity to the frozen discovery ranking, not an independent
# ligand-significance test. Shared prior remains a limitation.
comparison <- broad$comparison
if (!all(c("test_ligand", "discovery_broad_rank") %in% names(comparison))) {
  stop("Unexpected frozen broad comparison structure")
}
discovery_rank <- comparison$discovery_broad_rank[match(ligands, comparison$test_ligand)]
if (anyNA(discovery_rank)) stop("Frozen discovery ranking missing ligand(s)")
rank_rho <- vapply(seq_len(n_subsets), function(i) {
  suppressWarnings(cor(discovery_rank, ranks[i, ], method = "spearman"))
}, numeric(1))

subset_summary <- data.frame(
  subset_id = seq_len(n_subsets),
  significant_ligands_FDR05 = rowSums(fdr < 0.05),
  positive_corrected_AUPR_ligands = rowSums(subset_scores > 0),
  discovery_rank_spearman_rho = rank_rho,
  stringsAsFactors = FALSE
)

all_ligand_summary <- data.frame(
  test_ligand = ligands,
  median_aupr_corrected = apply(subset_scores, 2, median),
  median_rank = apply(ranks, 2, median),
  top5_fraction = colMeans(ranks <= 5),
  top10_fraction = colMeans(ranks <= 10),
  positive_aupr_fraction = colMeans(subset_scores > 0),
  empirical_FDR05_fraction = colMeans(fdr < 0.05),
  median_empirical_P = apply(pvals, 2, median),
  median_empirical_FDR = apply(fdr, 2, median),
  stringsAsFactors = FALSE
)

main <- c("SEMA7A", "EDIL3", "FLRT3", "PGF", "SPP1")
main_summary <- all_ligand_summary[match(main, all_ligand_summary$test_ligand), , drop = FALSE]
if (anyNA(main_summary$test_ligand)) stop("Main candidate missing from broad universe")

# Long-form per-subset results facilitate independent figure generation.
long <- data.frame(
  subset_id = rep(seq_len(n_subsets), each = length(ligands)),
  test_ligand = rep(ligands, times = n_subsets),
  aupr_corrected = as.vector(t(subset_scores)),
  rank = as.vector(t(ranks)),
  empirical_P = as.vector(t(pvals)),
  empirical_FDR = as.vector(t(fdr)),
  stringsAsFactors = FALSE
)

bin_audit <- data.frame(
  expression_decile = seq_len(n_bins),
  independent_receiver_697 = receiver_bin_n,
  size_matched_22_quota = quota,
  background_n = tabulate(bins[background], nbins = n_bins)
)

summary_stats <- data.frame(
  metric = c("n_subsets", "n_null", "receiver_full", "receiver_matched",
             "background", "ligands", "median_significant_ligands_per_subset",
             "mean_significant_ligands_per_subset", "fraction_subsets_with_any_FDR05",
             "fraction_subsets_with_zero_FDR05", "median_discovery_rank_spearman_rho"),
  value = c(n_subsets, n_null, length(receiver), n_match, length(background),
            length(ligands), median(subset_summary$significant_ligands_FDR05),
            mean(subset_summary$significant_ligands_FDR05),
            mean(subset_summary$significant_ligands_FDR05 > 0),
            mean(subset_summary$significant_ligands_FDR05 == 0),
            median(rank_rho)),
  stringsAsFactors = FALSE
)

write.csv(bin_audit, file.path(outdir, "18_expression_bin_audit.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(outdir, "18_summary_statistics.csv"), row.names = FALSE)
write.csv(subset_summary, file.path(outdir, "18_receiver_subset_summary.csv"), row.names = FALSE)
write.csv(all_ligand_summary, file.path(outdir, "18_all42_ligand_stability.csv"), row.names = FALSE)
write.csv(main_summary, file.path(outdir, "18_main5_candidate_stability.csv"), row.names = FALSE)
write.csv(long, file.path(outdir, "18_per_subset_ligand_results.csv"), row.names = FALSE)
saveRDS(list(
  settings = list(n_subsets = n_subsets, n_null = n_null, n_match = n_match,
                  n_cores = n_cores, seed = seed, run_tag = run_tag,
                  prior_md5 = observed_md5, quota = quota),
  ligand_order = ligands, subset_gene_sets = subset_sets,
  null_gene_sets = null_sets, subset_scores = subset_scores,
  null_scores = null_scores, pvals = pvals, fdr = fdr, ranks = ranks,
  subset_summary = subset_summary, all_ligand_summary = all_ligand_summary,
  main_summary = main_summary, summary_stats = summary_stats
), file.path(outdir, "18_complete_analysis.rds"))
capture.output(sessionInfo(), file = file.path(outdir, "18_sessionInfo.txt"))

cat("\n===== SIZE-MATCHED RESULTS =====\n")
print(summary_stats, row.names = FALSE)
cat("\n===== FIVE MAIN CANDIDATES =====\n")
print(main_summary, row.names = FALSE)
for (f in c(original_null_file, discovery_null_file)) {
  if (file.exists(f)) {
    x <- read.csv(f, check.names = FALSE)
    if ("empirical_FDR" %in% names(x)) {
      cat("Original frozen reference:", f, "=>", sum(x$empirical_FDR < 0.05, na.rm = TRUE),
          "/", nrow(x), "FDR < 0.05\n")
    }
  }
}
cat("\nOutputs saved in:", outdir, "\n")
cat("INTERPRETATION: This tests sensitivity to receiver size and sampling within\n")
cat("the independent 697-gene program; it does NOT validate direct signaling.\n")
cat("The 22-gene subsets are not independent biological replication cohorts.\n")
cat("SCRIPT 18 COMPLETED SUCCESSFULLY\n")
