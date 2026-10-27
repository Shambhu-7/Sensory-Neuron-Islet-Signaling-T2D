#!/usr/bin/env Rscript
# SCRIPT 20 — NicheNet same-expression-bin gene replacement sensitivity
# Non-destructive exploratory analysis following Scripts 18 and 19.
# For four post-hoc exceptional 22-gene receiver subsets, replace each
# high-prior-weight gene (and a lower-weight control) with genes from the
# SAME expression decile in the INDEPENDENT 697-gene receiver, excluding the
# original 22 genes. Every replacement retains exactly 22 genes and the
# original 10-bin quota. All 42 ligand activities use the EXACT scoring
# function and shared 22-gene expression-matched null from Script 18.
# Results are sensitivity diagnostics, NOT biological or causal validation.
# Run from project root. Does not alter frozen inputs or earlier results.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(nichenetr))

positive_int <- function(key, default) {
  x <- suppressWarnings(as.integer(Sys.getenv(key, as.character(default))))
  if (length(x) != 1L || is.na(x) || x < 1L) stop("Invalid ", key)
  x
}
n_replacements <- positive_int("N_REPLACEMENTS", 50L)
n_cores <- positive_int("N_CORES", 4L)
chunk_size <- positive_int("CHUNK_SIZE", 25L)
run_tag <- Sys.getenv("RUN_TAG", "pilot1")
if (!grepl("^[A-Za-z0-9_-]+$", run_tag)) stop("Unsafe RUN_TAG")
seed <- 200221156L

lt_file <- "data/raw/NicheNet/ligand_target_matrix_nsga2r_final.rds"
de_file <- "results/validation/GSE221156/Table_P3_GSE221156_Beta_T2D_vs_Normal_edgeR.csv"
independent_file <- "data/processed/GSE221156_Independent_NicheNet_16a.rds"
broad_file <- "data/processed/BroadLigandUniverse_NicheNet_16b.rds"
pilot_file <- "results/validation/GSE221156/NicheNet_receiver_size_matched/pilot1/18_complete_analysis.rds"
script19_file <- "results/validation/GSE221156/NicheNet_LOO_influence/pilot1/19_complete_analysis.rds"
outdir <- file.path("results/validation/GSE221156/NicheNet_same_bin_replacement", run_tag)
required <- c(lt_file, de_file, independent_file, broad_file, pilot_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required inputs:\n", paste(missing, collapse = "\n"))
if (file.exists(file.path(outdir, "20_complete_analysis.rds"))) {
  stop("Completed analysis already exists in ", outdir, "; choose a new RUN_TAG")
}

expected_md5 <- "b09606b04b2d4490418d9028c0e58b9f"
prior_md5 <- unname(tools::md5sum(lt_file))
if (!identical(prior_md5, expected_md5)) stop("Frozen NicheNet prior MD5 mismatch")
lt <- readRDS(lt_file)
if (!identical(dim(lt), c(33354L, 1226L))) stop("Unexpected prior dimensions")
independent <- readRDS(independent_file)
broad <- readRDS(broad_file)
pilot <- readRDS(pilot_file)
ligands <- sort(unique(as.character(broad$broad_ligands)))
if (length(ligands) != 42L || !identical(ligands, pilot$ligand_order) ||
    !all(ligands %in% colnames(lt))) stop("Frozen ligand universe mismatch")
receiver <- unique(as.character(independent$receiver_nichenet))
background <- unique(as.character(independent$background_nichenet))
if (length(receiver) != 697L || length(background) != 13712L ||
    !all(receiver %in% background) || !all(background %in% rownames(lt))) {
  stop("Frozen independent receiver/background mismatch")
}
if (!is.matrix(pilot$subset_scores) || !is.matrix(pilot$null_scores) ||
    !identical(colnames(pilot$subset_scores), ligands) ||
    !identical(colnames(pilot$null_scores), ligands) ||
    !identical(colnames(pilot$fdr), ligands) ||
    nrow(pilot$null_scores) != pilot$settings$n_null ||
    length(pilot$subset_gene_sets) < 100L ||
    !identical(as.integer(pilot$settings$n_match), 22L)) {
  stop("Script 18 pilot RDS structure mismatch")
}
if (any(!is.finite(pilot$null_scores)) || any(!is.finite(pilot$subset_scores))) {
  stop("Nonfinite Script 18 NicheNet scores")
}
ids <- c(21L, 40L, 94L, 100L)
expected_counts <- c(39L, 38L, 31L, 34L)
observed_counts <- as.integer(pilot$subset_summary$significant_ligands_FDR05[ids])
if (!identical(observed_counts, expected_counts)) {
  stop("Unexpected Script 18 exceptional subset counts")
}

de <- read.csv(de_file, check.names = FALSE)
if (!all(c("gene", "logCPM") %in% names(de))) stop("DE table missing gene/logCPM")
de$gene <- trimws(as.character(de$gene))
de <- de[!is.na(de$gene) & nzchar(de$gene) & !duplicated(de$gene), , drop = FALSE]
if (!setequal(intersect(de$gene, rownames(lt)), background)) {
  stop("DE background does not reproduce frozen background")
}
expr <- de$logCPM
names(expr) <- de$gene
expr <- expr[background]
if (any(!is.finite(expr))) stop("Nonfinite logCPM values")

# Exact expression-decile assignment used in Scripts 18 and 19.
make_expression_bins <- function(background_genes, expression_values, n_bins = 10L) {
  e <- expression_values[background_genes]
  if (anyNA(e)) stop("Missing expression")
  ord <- order(e, names(e))
  n <- length(ord)
  bin_ord <- pmin(n_bins, ceiling(seq_len(n) * n_bins / n))
  bins <- integer(n)
  bins[ord] <- bin_ord
  names(bins) <- names(e)
  bins
}
bins <- make_expression_bins(background, expr)
quota22 <- as.integer(pilot$settings$quota)
if (length(quota22) != 10L || sum(quota22) != 22L) stop("Invalid Script 18 quota")

# Restrict prior exactly as in Script 18. This is NOT a proxy score.
lt_sub <- lt[background, ligands, drop = FALSE]
rm(lt)
invisible(gc(verbose = FALSE))
run_activity <- function(geneset) {
  a <- nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = lt_sub,
    potential_ligands = ligands
  )
  v <- a$aupr_corrected
  names(v) <- as.character(a$test_ligand)
  v <- as.numeric(v[ligands])
  if (length(v) != length(ligands) || any(!is.finite(v))) {
    stop("Invalid NicheNet corrected AUPR")
  }
  v
}

# If Script 19 exists, confirm the exact influential genes it identified.
expected_top <- c(`21` = "TNFRSF10B", `40` = "CDKN1A",
                  `94` = "FAS", `100` = "ODC1")
if (file.exists(script19_file)) {
  s19 <- readRDS(script19_file)
  for (id in ids) {
    x <- s19$summary
    x <- x[x$subset_id == id & x$deletion_type == "top_weight", , drop = FALSE]
    if (nrow(x) != 1L || !identical(as.character(x$removed_gene),
                                   unname(expected_top[as.character(id)]))) {
      stop("Script 19 top gene mismatch for subset ", id)
    }
  }
}

# Manifest: one high-weight target and one lower-weight control per subset.
# For 21 and 94, controls are in the SAME bin as the target (as in Script 19).
# For 40 and 100, no same-bin internal control exists; choose the lowest-weight
# other gene in the 22-gene subset, label it explicitly as a DIFFERENT-bin
# control, and replace it with a gene from ITS OWN expression bin.
manifest <- list()
for (id in ids) {
  genes <- pilot$subset_gene_sets[[id]]
  if (length(genes) != 22L || anyDuplicated(genes) || !all(genes %in% receiver)) {
    stop("Invalid 22-gene subset ", id)
  }
  if (!identical(as.integer(tabulate(bins[genes], nbins = 10L)), quota22)) {
    stop("Script 18 expression quota mismatch for subset ", id)
  }
  weights <- rowMeans(lt_sub[genes, , drop = FALSE])
  names(weights) <- genes
  top_gene <- names(weights)[which.max(weights)]
  if (!identical(top_gene, unname(expected_top[as.character(id)]))) {
    stop("Highest-prior-weight gene mismatch for subset ", id, ": ", top_gene)
  }
  same_bin <- setdiff(genes[bins[genes] == bins[top_gene]], top_gene)
  if (length(same_bin)) {
    control_gene <- same_bin[which.min(weights[same_bin])]
    control_type <- "same_bin_low_weight_control"
  } else {
    others <- setdiff(genes, top_gene)
    control_gene <- others[which.min(weights[others])]
    control_type <- "other_bin_low_weight_control"
  }
  for (j in seq_along(c(top_gene, control_gene))) {
    removed <- c(top_gene, control_gene)[j]
    typ <- c("top_weight", control_type)[j]
    pool <- setdiff(receiver[bins[receiver] == bins[removed]], genes)
    pool <- sort(unique(pool))
    if (!length(pool)) stop("No replacement candidates for ", removed)
    # Sample without replacement within this target/control test.
    n_draw <- min(n_replacements, length(pool))
    manifest[[length(manifest) + 1L]] <- list(
      subset_id = id, deletion_type = typ, removed_gene = removed,
      removed_gene_bin = unname(bins[removed]),
      removed_gene_weight = unname(weights[removed]),
      removed_gene_weight_share = unname(weights[removed] / sum(weights)),
      original_genes = genes, available_pool = pool, n_draw = n_draw
    )
  }
}

# Deterministic manifest, independent of chunk size and worker count.
set.seed(seed)
for (i in seq_along(manifest)) {
  manifest[[i]]$replacement_genes <- sample(manifest[[i]]$available_pool,
                                            manifest[[i]]$n_draw, replace = FALSE)
}
manifest_csv <- do.call(rbind, lapply(seq_along(manifest), function(i) {
  z <- manifest[[i]]
  data.frame(test_id = i, subset_id = z$subset_id,
             deletion_type = z$deletion_type, removed_gene = z$removed_gene,
             removed_gene_bin = z$removed_gene_bin,
             removed_gene_weight_share = z$removed_gene_weight_share,
             n_available_replacements = length(z$available_pool),
             n_sampled_replacements = z$n_draw, stringsAsFactors = FALSE)
}))

jobs <- list()
for (i in seq_along(manifest)) {
  z <- manifest[[i]]
  for (k in seq_along(z$replacement_genes)) {
    replacement <- z$replacement_genes[k]
    genes22 <- c(z$original_genes[z$original_genes != z$removed_gene], replacement)
    if (length(genes22) != 22L || anyDuplicated(genes22) ||
        !all(genes22 %in% receiver) ||
        !identical(as.integer(tabulate(bins[genes22], nbins = 10L)), quota22)) {
      stop("Replacement does not preserve receiver size and expression quotas")
    }
    jobs[[length(jobs) + 1L]] <- list(
      test_id = i, subset_id = z$subset_id, deletion_type = z$deletion_type,
      removed_gene = z$removed_gene, replacement_gene = replacement,
      replacement_bin = unname(bins[replacement]), genes22 = genes22
    )
  }
}

# Include exact source hashes, package version, and the full deterministic job
# manifest in checkpoint identity. Never reuse an incompatible checkpoint.
input_hashes <- unname(tools::md5sum(required))
job_string <- paste(vapply(jobs, function(z) {
  paste(z$test_id, z$removed_gene, z$replacement_gene, sep = ":")
}, character(1)), collapse = "|")
fingerprint <- paste(c(input_hashes, as.character(packageVersion("nichenetr")),
                       seed, n_replacements, paste(ligands, collapse = ","),
                       job_string), collapse = "|")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path(outdir, "checkpoint_replacement_scores.rds")
scores <- matrix(numeric(), nrow = 0L, ncol = length(ligands),
                 dimnames = list(NULL, ligands))
if (file.exists(checkpoint)) {
  old <- readRDS(checkpoint)
  if (!identical(old$fingerprint, fingerprint)) {
    stop("Checkpoint mismatch: choose a new RUN_TAG; do not delete old results")
  }
  scores <- old$scores
  if (!is.matrix(scores) || ncol(scores) != length(ligands) ||
      !identical(colnames(scores), ligands) || nrow(scores) > length(jobs) ||
      any(!is.finite(scores))) stop("Malformed replacement checkpoint")
}

cat("\nSCRIPT 20 — SAME-BIN RECEIVER-GENE REPLACEMENT SENSITIVITY\n")
cat("Frozen prior MD5:", prior_md5, "\n")
cat("Frozen Script 18 null sets:", nrow(pilot$null_scores), "(REUSED, not recomputed)\n")
cat("Original receiver size: 22 | expression quota:", paste(quota22, collapse = ","), "\n")
cat("Exceptional subsets:", paste(ids, collapse = ","), "\n")
cat("Tests:", length(manifest), "| replacement jobs:", length(jobs),
    "| cores:", n_cores, "| chunk size:", chunk_size, "\n")
cat("Output:", outdir, "\n")
print(manifest_csv, row.names = FALSE)
if (nrow(scores)) cat("Resuming from", nrow(scores), "/", length(jobs), "scored replacements\n")

if (nrow(scores) < length(jobs)) {
  for (start in seq.int(nrow(scores) + 1L, length(jobs), by = chunk_size)) {
    idx <- seq.int(start, min(start + chunk_size - 1L, length(jobs)))
    t0 <- proc.time()[["elapsed"]]
    ans <- parallel::mclapply(idx, function(i) run_activity(jobs[[i]]$genes22),
                              mc.cores = n_cores, mc.preschedule = TRUE,
                              mc.set.seed = FALSE)
    bad <- vapply(ans, function(z) inherits(z, "try-error") ||
                    !is.numeric(z) || length(z) != length(ligands) ||
                    any(!is.finite(z)), logical(1))
    if (any(bad)) stop("NicheNet worker failed at replacement job ", start)
    block <- do.call(rbind, ans)
    colnames(block) <- ligands
    scores <- rbind(scores, block)
    saveRDS(list(fingerprint = fingerprint, scores = scores), checkpoint)
    cat(sprintf("replacements: %d/%d | batch %.1f sec\n", nrow(scores),
                length(jobs), proc.time()[["elapsed"]] - t0))
  }
}

# Compare all replacement scores against the SAME 22-gene null used in
# Script 18. Empirical one-sided P and BH across all 42 ligands per set.
null <- pilot$null_scores[, ligands, drop = FALSE]
pvals <- matrix(NA_real_, nrow = length(jobs), ncol = length(ligands),
                dimnames = list(NULL, ligands))
fdr <- pvals
for (i in seq_along(jobs)) {
  exceed <- vapply(seq_along(ligands), function(j) {
    sum(null[, j] >= scores[i, j])
  }, integer(1))
  pvals[i, ] <- (1 + exceed) / (1 + nrow(null))
  fdr[i, ] <- p.adjust(pvals[i, ], method = "BH")
}

main5 <- c("SEMA7A", "EDIL3", "FLRT3", "PGF", "SPP1")
if (!all(main5 %in% ligands)) stop("Main-five ligands missing")
rows <- vector("list", length(jobs))
job_summary <- vector("list", length(jobs))
for (i in seq_along(jobs)) {
  z <- jobs[[i]]
  orig <- as.numeric(pilot$subset_scores[z$subset_id, ligands])
  old_fdr <- as.numeric(pilot$fdr[z$subset_id, ligands])
  new <- as.numeric(scores[i, ])
  new_fdr <- as.numeric(fdr[i, ])
  rank_orig <- rank(-orig, ties.method = "average")
  rank_new <- rank(-new, ties.method = "average")
  rho <- suppressWarnings(cor(rank_orig, rank_new, method = "spearman"))
  rows[[i]] <- data.frame(
    job_id = i, test_id = z$test_id, subset_id = z$subset_id,
    deletion_type = z$deletion_type, removed_gene = z$removed_gene,
    replacement_gene = z$replacement_gene, expression_bin = z$replacement_bin,
    test_ligand = ligands, original_22_AUPR = orig, replacement_22_AUPR = new,
    delta_AUPR = new - orig, original_22_FDR = old_fdr,
    replacement_22_P = as.numeric(pvals[i, ]), replacement_22_FDR = new_fdr,
    original_significant = old_fdr < 0.05,
    replacement_significant = new_fdr < 0.05,
    original_rank = rank_orig, replacement_rank = rank_new,
    stringsAsFactors = FALSE
  )
  job_summary[[i]] <- data.frame(
    job_id = i, test_id = z$test_id, subset_id = z$subset_id,
    deletion_type = z$deletion_type, removed_gene = z$removed_gene,
    replacement_gene = z$replacement_gene, expression_bin = z$replacement_bin,
    significant_original = sum(old_fdr < 0.05),
    significant_replacement = sum(new_fdr < 0.05),
    delta_significant = sum(new_fdr < 0.05) - sum(old_fdr < 0.05),
    median_delta_AUPR = median(new - orig),
    mean_abs_delta_AUPR = mean(abs(new - orig)),
    spearman_rank = rho,
    main5_significant_original = sum(old_fdr[match(main5, ligands)] < 0.05),
    main5_significant_replacement = sum(new_fdr[match(main5, ligands)] < 0.05),
    stringsAsFactors = FALSE
  )
}
all42 <- do.call(rbind, rows)
job_table <- do.call(rbind, job_summary)

# One row per tested gene. Descriptive fractions, not permutation P-values:
# targets were chosen after inspecting the original exceptional subsets.
target_summary <- do.call(rbind, lapply(seq_along(manifest), function(i) {
  z <- manifest[[i]]
  x <- job_table[job_table$test_id == i, , drop = FALSE]
  data.frame(
    test_id = i, subset_id = z$subset_id, deletion_type = z$deletion_type,
    removed_gene = z$removed_gene,
    removed_gene_expression_bin = z$removed_gene_bin,
    removed_gene_prior_weight_share = z$removed_gene_weight_share,
    n_replacements = nrow(x), significant_original = x$significant_original[1L],
    median_significant_replacement = median(x$significant_replacement),
    range_significant_min = min(x$significant_replacement),
    range_significant_max = max(x$significant_replacement),
    fraction_zero_significant = mean(x$significant_replacement == 0L),
    fraction_any_significant = mean(x$significant_replacement > 0L),
    median_delta_significant = median(x$delta_significant),
    median_delta_AUPR = median(x$median_delta_AUPR),
    median_rank_spearman = median(x$spearman_rank, na.rm = TRUE),
    median_main5_significant_replacement = median(x$main5_significant_replacement),
    stringsAsFactors = FALSE
  )
}))

# Ligand-level retention frequency across replacements, for all 42 ligands.
ligand_stability <- do.call(rbind, lapply(seq_along(manifest), function(i) {
  x <- all42[all42$test_id == i, , drop = FALSE]
  z <- manifest[[i]]
  data.frame(
    test_id = i, subset_id = z$subset_id, deletion_type = z$deletion_type,
    removed_gene = z$removed_gene, test_ligand = ligands,
    original_22_FDR = as.numeric(pilot$fdr[z$subset_id, ligands]),
    fraction_significant_after_replacement = vapply(ligands, function(g) {
      mean(x$replacement_significant[x$test_ligand == g])
    }, numeric(1)),
    median_delta_AUPR = vapply(ligands, function(g) {
      median(x$delta_AUPR[x$test_ligand == g])
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}))

write.csv(manifest_csv, file.path(outdir, "20_test_manifest.csv"), row.names = FALSE)
write.csv(job_table, file.path(outdir, "20_replacement_level_summary.csv"), row.names = FALSE)
write.csv(target_summary, file.path(outdir, "20_target_summary.csv"), row.names = FALSE)
write.csv(all42, file.path(outdir, "20_all42_ligand_replacements.csv"), row.names = FALSE)
write.csv(ligand_stability, file.path(outdir, "20_all42_ligand_stability.csv"), row.names = FALSE)
write.csv(ligand_stability[ligand_stability$test_ligand %in% main5, , drop = FALSE],
          file.path(outdir, "20_main5_ligand_stability.csv"), row.names = FALSE)
write.csv(data.frame(expression_decile = seq_len(10L),
                     original_22_quota = quota22,
                     receiver_697_bin_n = tabulate(bins[receiver], nbins = 10L),
                     background_bin_n = tabulate(bins[background], nbins = 10L)),
          file.path(outdir, "20_expression_bin_audit.csv"), row.names = FALSE)
saveRDS(list(
  settings = list(n_replacements = n_replacements, n_cores = n_cores,
                  chunk_size = chunk_size, seed = seed, run_tag = run_tag,
                  prior_md5 = prior_md5, source_md5 = input_hashes,
                  original_null_n = nrow(null), original_quota = quota22),
  manifest = manifest, jobs = jobs, replacement_scores = scores,
  empirical_p = pvals, empirical_fdr = fdr, target_summary = target_summary,
  job_summary = job_table, ligand_stability = ligand_stability,
  original_pilot_path = pilot_file
), file.path(outdir, "20_complete_analysis.rds"))
capture.output(sessionInfo(), file = file.path(outdir, "20_sessionInfo.txt"))

cat("\n===== SCRIPT 20 SUMMARY =====\n")
print(target_summary[, c("subset_id", "deletion_type", "removed_gene",
                         "n_replacements", "significant_original",
                         "median_significant_replacement", "range_significant_min",
                         "range_significant_max", "fraction_zero_significant",
                         "median_delta_AUPR", "median_rank_spearman")],
      row.names = FALSE)
cat("\nOutputs saved in:", outdir, "\n")
cat("INTERPRETATION: Descriptive, post-selection sensitivity analysis only.\n")
cat("22-gene receiver size and expression quota held constant for every replacement.\n")
cat("Original Script 18 null reused; no new null permutations required.\n")
cat("Replacement candidates come from the same independent 697-gene receiver.\n")
cat("No inference of direct signaling, causality, or cohort-level replication.\n")
cat("SCRIPT 20 COMPLETED SUCCESSFULLY\n")
