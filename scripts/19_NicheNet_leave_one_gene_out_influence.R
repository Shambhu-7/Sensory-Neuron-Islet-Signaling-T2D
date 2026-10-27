#!/usr/bin/env Rscript
# SCRIPT 19 — NicheNet leave-one-gene-out receiver influence analysis
# Non-destructive sensitivity analysis following Script 18.
# Targeted mode: delete highest-prior-weight gene in each exceptional 22-gene
# receiver subset, and an expression-decile-matched low-weight control gene.
# Every 21-gene test uses a NEW 21-gene expression-matched empirical null.
# Null sets are shared only when the exact 10-bin quota is identical.
# Exploratory, not evidence of causal ligand-receptor signaling.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(nichenetr))

positive_int <- function(key, default) {
  z <- suppressWarnings(as.integer(Sys.getenv(key, as.character(default))))
  if (length(z) != 1L || is.na(z) || z < 1L) stop("Invalid ", key)
  z
}
n_null <- positive_int("N_NULL", 1000L)
n_cores <- positive_int("N_CORES", 4L)
chunk_size <- positive_int("CHUNK_SIZE", 50L)
run_tag <- Sys.getenv("RUN_TAG", "pilot1")
if (!grepl("^[A-Za-z0-9_-]+$", run_tag)) stop("Unsafe RUN_TAG")
mode <- Sys.getenv("MODE", "targeted")
if (!mode %in% c("targeted", "all")) stop("MODE must be targeted or all")
if (n_null < 999L) warning("Coarse empirical P resolution; pilot only")

lt_file <- "data/raw/NicheNet/ligand_target_matrix_nsga2r_final.rds"
de_file <- "results/validation/GSE221156/Table_P3_GSE221156_Beta_T2D_vs_Normal_edgeR.csv"
independent_file <- "data/processed/GSE221156_Independent_NicheNet_16a.rds"
broad_file <- "data/processed/BroadLigandUniverse_NicheNet_16b.rds"
pilot_file <- "results/validation/GSE221156/NicheNet_receiver_size_matched/pilot1/18_complete_analysis.rds"
required <- c(lt_file, de_file, independent_file, broad_file, pilot_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing inputs:\n", paste(missing, collapse = "\n"))

expected_md5 <- "b09606b04b2d4490418d9028c0e58b9f"
actual_md5 <- unname(tools::md5sum(lt_file))
if (!identical(actual_md5, expected_md5)) stop("Frozen prior MD5 mismatch")
lt <- readRDS(lt_file)
if (!identical(dim(lt), c(33354L, 1226L))) stop("Unexpected prior dimensions")
independent <- readRDS(independent_file)
broad <- readRDS(broad_file)
pilot <- readRDS(pilot_file)
ligands <- sort(unique(as.character(broad$broad_ligands)))
if (length(ligands) != 42L || !identical(ligands, pilot$ligand_order) ||
    !all(ligands %in% colnames(lt))) stop("Ligand universe mismatch")
receiver <- unique(as.character(independent$receiver_nichenet))
background <- unique(as.character(independent$background_nichenet))
if (length(receiver) != 697L || length(background) != 13712L ||
    !all(receiver %in% background) || !all(background %in% rownames(lt))) {
  stop("Frozen independent receiver/background mismatch")
}
if (ncol(pilot$subset_scores) != 42L ||
    !identical(colnames(pilot$subset_scores), ligands) ||
    !identical(colnames(pilot$fdr), ligands) ||
    length(pilot$subset_gene_sets) < 100L ||
    !all(vapply(pilot$subset_gene_sets, length, integer(1)) == 22L)) {
  stop("Script 18 pilot RDS has unexpected structure")
}

de <- read.csv(de_file, check.names = FALSE)
if (!all(c("gene", "logCPM") %in% names(de))) stop("DE input columns missing")
de$gene <- trimws(as.character(de$gene))
de <- de[!is.na(de$gene) & nzchar(de$gene) & !duplicated(de$gene), , drop = FALSE]
if (!setequal(intersect(de$gene, rownames(lt)), background)) {
  stop("Expression background does not match frozen input")
}
expr <- de$logCPM
names(expr) <- de$gene
expr <- expr[background]
if (any(!is.finite(expr))) stop("Nonfinite background expression")

# Exactly reproduce Script 18's expression-decile algorithm.
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
quota22 <- tabulate(bins[pilot$subset_gene_sets[[21]]], nbins = 10L)
if (!identical(as.integer(quota22), as.integer(pilot$settings$quota))) {
  stop("Script 18 expression bin quotas do not reproduce")
}

lt_sub <- lt[background, ligands, drop = FALSE]
rm(lt)
gc(verbose = FALSE)
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
  if (length(v) != 42L || any(!is.finite(v))) stop("Invalid NicheNet activity")
  v
}

# The 4 exceptional subsets are defined by the frozen Script 18 pilot.
ids <- c(21L, 40L, 94L, 100L)
if (!identical(as.integer(pilot$subset_summary$significant_ligands_FDR05[ids]),
               c(39L, 38L, 31L, 34L))) {
  stop("Script 18 exceptional subset counts do not match audited pilot")
}

# Build an auditable test manifest. In targeted mode, choose a same-bin control
# with the smallest mean prior weight. If no alternative exists, skip that
# control rather than silently compare against a different bin.
manifest <- list()
for (id in ids) {
  genes <- pilot$subset_gene_sets[[id]]
  if (length(genes) != 22L || anyDuplicated(genes) || !all(genes %in% receiver)) {
    stop("Invalid receiver subset ", id)
  }
  w <- rowMeans(lt_sub[genes, , drop = FALSE])
  names(w) <- genes
  top_gene <- names(w)[which.max(w)]
  if (mode == "all") {
    choices <- genes
    types <- ifelse(genes == top_gene, "top_weight", "other_gene")
  } else {
    same_bin <- genes[bins[genes] == bins[top_gene]]
    control_candidates <- setdiff(same_bin, top_gene)
    if (!length(control_candidates)) {
      warning("No same-expression-bin control for subset ", id,
              "; testing top gene only")
      choices <- top_gene
      types <- "top_weight"
    } else {
      control_gene <- control_candidates[which.min(w[control_candidates])]
      choices <- c(top_gene, control_gene)
      types <- c("top_weight", "same_bin_low_weight_control")
    }
  }
  for (k in seq_along(choices)) {
    g <- choices[k]
    genes21 <- setdiff(genes, g)
    q21 <- tabulate(bins[genes21], nbins = 10L)
    if (length(genes21) != 21L || sum(q21) != 21L) stop("Invalid 21-gene test")
    manifest[[length(manifest) + 1L]] <- list(
      subset_id = id, removed_gene = g, deletion_type = types[k],
      removed_gene_prior_weight = unname(w[g]),
      removed_gene_weight_share = unname(w[g] / sum(w)),
      removed_gene_expression_bin = unname(bins[g]),
      genes21 = genes21, quota21 = q21,
      quota_key = paste(q21, collapse = "_")
    )
  }
}

outdir <- file.path("results/validation/GSE221156/NicheNet_LOO_influence", run_tag)
if (dir.exists(outdir) && file.exists(file.path(outdir, "19_complete_analysis.rds"))) {
  stop("Output directory already contains a completed analysis. Choose a new RUN_TAG")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sample_stratified <- function(pool, bins, quota) {
  z <- unlist(lapply(seq_along(quota), function(b) {
    if (!quota[b]) return(character())
    available <- pool[bins[pool] == b]
    if (length(available) < quota[b]) stop("Insufficient genes in bin ", b)
    sample(available, quota[b], replace = FALSE)
  }), use.names = FALSE)
  if (length(z) != sum(quota) || anyDuplicated(z)) stop("Invalid sampled null")
  z
}

# The null is cached by exact 10-bin quota. Every test sharing a quota
# has exactly the same null distribution; all tests have 21 genes.
keys <- unique(vapply(manifest, `[[`, character(1), "quota_key"))
null_by_key <- list()
base_fingerprint <- paste(
  unname(tools::md5sum(required)),
  as.character(packageVersion("nichenetr")),
  n_null, paste(ligands, collapse = ","), sep = "|"
)

cat("\nSCRIPT 19 — LEAVE-ONE-GENE-OUT NICHE-NET INFLUENCE\n")
cat("Mode:", mode, "| tests:", length(manifest),
    "| distinct 21-gene null quotas:", length(keys),
    "| nulls per quota:", n_null, "| cores:", n_cores, "\n")
cat("Frozen prior MD5:", actual_md5, "\n")
cat("Output:", outdir, "\n")

for (key_index in seq_along(keys)) {
  key <- keys[key_index]
  q <- manifest[[which(vapply(manifest, function(z) identical(z$quota_key, key), logical(1)))[1L]]]$quota21
  fingerprint <- paste(base_fingerprint, key, key_index, sep = "|")
  checkpoint <- file.path(outdir, paste0("checkpoint_null_quota_", key, ".rds"))
  # Use deterministic seed per quota index (independent of cores/chunks).
  set.seed(190221156L + key_index)
  null_sets <- replicate(n_null, sample_stratified(background, bins, q), simplify = FALSE)
  scores <- matrix(numeric(), nrow = 0L, ncol = length(ligands),
                   dimnames = list(NULL, ligands))
  if (file.exists(checkpoint)) {
    old <- readRDS(checkpoint)
    if (!identical(old$fingerprint, fingerprint)) {
      stop("Checkpoint settings mismatch: ", checkpoint, "; use new RUN_TAG")
    }
    scores <- old$scores
    if (!is.matrix(scores) || ncol(scores) != 42L ||
        !identical(colnames(scores), ligands) || nrow(scores) > n_null ||
        any(!is.finite(scores))) stop("Malformed checkpoint: ", checkpoint)
    cat("Resuming quota ", key_index, " from ", nrow(scores), "/", n_null, "\n", sep = "")
  }
  if (nrow(scores) < n_null) {
    for (start in seq.int(nrow(scores) + 1L, n_null, by = chunk_size)) {
      idx <- seq.int(start, min(start + chunk_size - 1L, n_null))
      t0 <- proc.time()[["elapsed"]]
      ans <- parallel::mclapply(idx, function(i) run_activity(null_sets[[i]]),
                                mc.cores = n_cores, mc.preschedule = TRUE,
                                mc.set.seed = FALSE)
      bad <- vapply(ans, function(z) inherits(z, "try-error") ||
                      !is.numeric(z) || length(z) != length(ligands) ||
                      any(!is.finite(z)), logical(1))
      if (any(bad)) stop("NicheNet worker failure in quota ", key, " chunk ", start)
      chunk <- do.call(rbind, ans)
      colnames(chunk) <- ligands
      scores <- rbind(scores, chunk)
      saveRDS(list(fingerprint = fingerprint, scores = scores), checkpoint)
      cat(sprintf("null quota %d/%d: %d/%d | %.1f sec\n",
                  key_index, length(keys), nrow(scores), n_null,
                  proc.time()[["elapsed"]] - t0))
    }
  }
  null_by_key[[key]] <- scores
}

# Observed activity for each distinct 21-gene receiver set.
observed21 <- parallel::mclapply(seq_along(manifest),
  function(i) run_activity(manifest[[i]]$genes21),
  mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
if (any(vapply(observed21, function(z) inherits(z, "try-error") ||
               !is.numeric(z) || length(z) != 42L || any(!is.finite(z)), logical(1)))) {
  stop("NicheNet observed receiver scoring failed")
}

rows <- vector("list", length(manifest))
summary_rows <- vector("list", length(manifest))
main <- c("SEMA7A", "EDIL3", "FLRT3", "PGF", "SPP1")
for (i in seq_along(manifest)) {
  z <- manifest[[i]]
  id <- z$subset_id
  obs <- observed21[[i]]
  null <- null_by_key[[z$quota_key]]
  exceed <- vapply(seq_along(ligands), function(j) sum(null[, j] >= obs[j]), integer(1))
  p21 <- (1 + exceed) / (1 + nrow(null))
  fdr21 <- p.adjust(p21, method = "BH")
  old_scores <- as.numeric(pilot$subset_scores[id, ligands])
  old_fdr <- as.numeric(pilot$fdr[id, ligands])
  rank22 <- rank(-old_scores, ties.method = "average")
  rank21 <- rank(-obs, ties.method = "average")
  rho <- suppressWarnings(cor(rank22, rank21, method = "spearman"))
  rows[[i]] <- data.frame(
    subset_id = id, deletion_type = z$deletion_type,
    removed_gene = z$removed_gene,
    removed_gene_expression_bin = z$removed_gene_expression_bin,
    removed_gene_prior_weight = z$removed_gene_prior_weight,
    removed_gene_weight_share = z$removed_gene_weight_share,
    test_ligand = ligands,
    aupr_22 = old_scores, aupr_21 = obs,
    delta_aupr_21_minus_22 = obs - old_scores,
    rank_22 = rank22, rank_21 = rank21,
    empirical_FDR_22 = old_fdr,
    empirical_P_21 = p21, empirical_FDR_21 = fdr21,
    significant_22 = old_fdr < 0.05,
    significant_21 = fdr21 < 0.05,
    quota_21 = z$quota_key,
    stringsAsFactors = FALSE
  )
  summary_rows[[i]] <- data.frame(
    subset_id = id, deletion_type = z$deletion_type,
    removed_gene = z$removed_gene,
    removed_gene_expression_bin = z$removed_gene_expression_bin,
    removed_gene_prior_weight = z$removed_gene_prior_weight,
    removed_gene_weight_share = z$removed_gene_weight_share,
    significant_ligands_22 = sum(old_fdr < 0.05),
    significant_ligands_21 = sum(fdr21 < 0.05),
    change_significant_21_minus_22 = sum(fdr21 < 0.05) - sum(old_fdr < 0.05),
    median_delta_AUPR = median(obs - old_scores),
    mean_abs_delta_AUPR = mean(abs(obs - old_scores)),
    spearman_rank_22_vs_21 = rho,
    n_main5_significant_22 = sum(old_fdr[match(main, ligands)] < 0.05),
    n_main5_significant_21 = sum(fdr21[match(main, ligands)] < 0.05),
    quota_21 = z$quota_key,
    stringsAsFactors = FALSE
  )
}
long <- do.call(rbind, rows)
summary <- do.call(rbind, summary_rows)
main5 <- long[long$test_ligand %in% main, , drop = FALSE]

# Write only within new Script 19 output directory.
write.csv(summary, file.path(outdir, "19_summary.csv"), row.names = FALSE)
write.csv(long, file.path(outdir, "19_all42_ligand_comparison.csv"), row.names = FALSE)
write.csv(main5, file.path(outdir, "19_main5_ligand_comparison.csv"), row.names = FALSE)
write.csv(data.frame(
  quota_key = keys,
  n_null = rep(n_null, length(keys)),
  stringsAsFactors = FALSE
), file.path(outdir, "19_null_quota_audit.csv"), row.names = FALSE)
saveRDS(list(settings = list(n_null = n_null, n_cores = n_cores, mode = mode,
                             run_tag = run_tag, prior_md5 = actual_md5),
             manifest = manifest, summary = summary, long = long,
             null_scores_by_quota = null_by_key,
             observed_21 = observed21,
             original_pilot_path = pilot_file),
        file.path(outdir, "19_complete_analysis.rds"))
capture.output(sessionInfo(), file = file.path(outdir, "19_sessionInfo.txt"))

cat("\n===== SCRIPT 19 SUMMARY =====\n")
print(summary[, c("subset_id", "deletion_type", "removed_gene",
                  "removed_gene_weight_share", "significant_ligands_22",
                  "significant_ligands_21", "change_significant_21_minus_22",
                  "median_delta_AUPR", "spearman_rank_22_vs_21")], row.names = FALSE)
cat("\nOutputs saved in:", outdir, "\n")
cat("CAVEAT: This is a receiver-gene influence sensitivity analysis.\n")
cat("It does not establish direct ligand-receptor signaling or causality.\n")
cat("SCRIPT 19 COMPLETED SUCCESSFULLY\n")
