#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

cat(strrep("=", 78), "\n", sep = "")
cat("GSE221156 — INDEPENDENT DONOR-LEVEL BETA-CELL VALIDATION\n")
cat(strrep("=", 78), "\n", sep = "")

pb_file <- "data/processed/GSE221156_Beta_DonorPseudobulk_RawCounts.tsv"
meta_file <- "results/validation/GSE221156/GSE221156_Donor_Metadata_Audit.csv"
receiver_file <- "results/islet/Table_P3_Beta_T2D_NicheNet_PrimaryReceiverGenes.csv"
core_file <- "results/islet/Table_P3_Beta_T2D_NicheNet_HighConfidenceCore.csv"
outdir <- "results/validation/GSE221156"
reprodir <- "reproducibility"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(reprodir, recursive = TRUE, showWarnings = FALSE)

required <- c(pb_file, meta_file, receiver_file, core_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required file(s):\n", paste(missing, collapse = "\n"))

cat("\nINPUT FILES: PASS\n")
for (f in required) cat("  ", f, "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. Read donor pseudobulk and collapse duplicate gene symbols by summing counts
# -----------------------------------------------------------------------------
pb <- read.delim(pb_file, check.names = FALSE, quote = "", comment.char = "")
if (!all(c("ensembl_id", "gene_symbol") %in% names(pb))) {
  stop("Pseudobulk table lacks ensembl_id/gene_symbol columns.")
}

donor_cols <- setdiff(names(pb), c("ensembl_id", "gene_symbol"))
counts0 <- as.matrix(pb[, donor_cols, drop = FALSE])
storage.mode(counts0) <- "numeric"
if (any(!is.finite(counts0)) || any(counts0 < 0)) stop("Invalid pseudobulk counts.")
if (!all(abs(counts0 - round(counts0)) < 1e-8)) stop("Pseudobulk counts are not integer-like.")

sym <- trimws(pb$gene_symbol)
ens <- trimws(pb$ensembl_id)
sym[is.na(sym) | sym == ""] <- ens[is.na(sym) | sym == ""]

# Sum the small number of duplicated symbols.
counts <- rowsum(counts0, group = sym, reorder = FALSE)
counts <- round(counts)
cat(sprintf("\nPseudobulk input: %s feature rows x %s donors\n", nrow(counts0), ncol(counts0)))
cat(sprintf("Unique gene symbols after duplicate collapse: %s\n", nrow(counts)))
cat(sprintf("Duplicate feature rows collapsed: %s\n", nrow(counts0) - nrow(counts)))

# -----------------------------------------------------------------------------
# 2. Donor metadata and prespecified primary comparison
# -----------------------------------------------------------------------------
meta <- read.csv(meta_file, check.names = FALSE)
if (!all(c("donor_id", "disease", "n_beta_cells") %in% names(meta))) {
  stop("Donor metadata lacks donor_id/disease/n_beta_cells.")
}
if (anyDuplicated(meta$donor_id)) stop("Duplicate donor IDs in donor metadata.")
if (!setequal(colnames(counts), meta$donor_id)) stop("Pseudobulk donor IDs do not match metadata.")
meta <- meta[match(colnames(counts), meta$donor_id), , drop = FALSE]
stopifnot(identical(colnames(counts), meta$donor_id))

cat("\nDONORS BY DISEASE (all 48)\n")
print(table(meta$disease))

keep_primary <- meta$disease %in% c("normal", "type 2 diabetes mellitus")
meta_p <- meta[keep_primary, , drop = FALSE]
counts_p <- counts[, keep_primary, drop = FALSE]
meta_p$group <- factor(ifelse(meta_p$disease == "type 2 diabetes mellitus", "T2D", "Normal"),
                       levels = c("Normal", "T2D"))

cat("\nPRIMARY COMPARISON: T2D vs normal\n")
cat(sprintf("  Normal donors: %d\n", sum(meta_p$group == "Normal")))
cat(sprintf("  T2D donors: %d\n", sum(meta_p$group == "T2D")))
cat(sprintf("  Total donors: %d\n", nrow(meta_p)))
cat(sprintf("  Beta cells represented: %s\n", format(sum(meta_p$n_beta_cells), big.mark = ",")))

# Metadata cross-tabs are descriptive and do not alter the prespecified primary model.
for (v in c("sex", "Center", "self_reported_ethnicity")) {
  if (v %in% names(meta_p)) {
    cat("\n", toupper(v), " x DISEASE\n", sep = "")
    print(table(meta_p[[v]], meta_p$group, useNA = "ifany"))
  }
}

# -----------------------------------------------------------------------------
# 3. Primary edgeR donor-level QL model: ~ group
# -----------------------------------------------------------------------------
y <- DGEList(counts = counts_p, samples = data.frame(meta_p, row.names = meta_p$donor_id))
keep_gene <- filterByExpr(y, group = meta_p$group)
y <- y[keep_gene, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ group, data = meta_p)
if (qr(design)$rank != ncol(design)) stop("Primary design matrix is not full rank.")

y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "groupT2D")
tab <- topTags(qlf, n = Inf, sort.by = "none")$table
tab$gene <- rownames(tab)
tab <- tab[, c("gene", setdiff(names(tab), "gene"))]
tab$direction <- ifelse(tab$logFC > 0, "T2D_up", ifelse(tab$logFC < 0, "T2D_down", "zero"))
tab$FDR05 <- tab$FDR < 0.05

de_file <- file.path(outdir, "Table_P3_GSE221156_Beta_T2D_vs_Normal_edgeR.csv")
write.csv(tab, de_file, row.names = FALSE)

cat("\nPRIMARY edgeR RESULTS\n")
cat(sprintf("  Genes before filter: %d\n", nrow(counts_p)))
cat(sprintf("  Genes retained: %d\n", nrow(y)))
cat(sprintf("  FDR < 0.05: %d\n", sum(tab$FDR < 0.05)))
cat(sprintf("  T2D-up at FDR < 0.05: %d\n", sum(tab$FDR < 0.05 & tab$logFC > 0)))
cat(sprintf("  T2D-down at FDR < 0.05: %d\n", sum(tab$FDR < 0.05 & tab$logFC < 0)))

# -----------------------------------------------------------------------------
# 4. Prespecified receiver-program validation
# -----------------------------------------------------------------------------
receiver <- read.csv(receiver_file, check.names = FALSE)
core <- read.csv(core_file, check.names = FALSE)
if (!all(c("gene", "primary_logFC") %in% names(receiver))) stop("Receiver table missing gene/primary_logFC.")
if (!"gene" %in% names(core)) stop("Core table missing gene.")

receiver <- receiver[!duplicated(receiver$gene), , drop = FALSE]
core_genes <- unique(core$gene)
receiver$original_direction <- ifelse(receiver$primary_logFC > 0, "T2D_up", "T2D_down")

val <- merge(receiver[, c("gene", "primary_logFC", "original_direction")],
             tab[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")],
             by = "gene", all.x = TRUE, sort = FALSE)
val$represented_in_GSE221156 <- val$gene %in% rownames(counts)
val$tested_by_edgeR <- !is.na(val$logFC)
val$validation_direction <- ifelse(is.na(val$logFC), NA,
                                   ifelse(val$logFC > 0, "T2D_up", "T2D_down"))
val$same_direction <- ifelse(val$tested_by_edgeR,
                             sign(val$primary_logFC) == sign(val$logFC), NA)
val$validation_P05 <- ifelse(val$tested_by_edgeR, val$PValue < 0.05, NA)
val$validation_FDR05 <- ifelse(val$tested_by_edgeR, val$FDR < 0.05, NA)
val$is_high_confidence_core <- val$gene %in% core_genes

receiver_out <- file.path(outdir, "Table_P3_GSE221156_27Gene_Receiver_Validation.csv")
write.csv(val, receiver_out, row.names = FALSE)

summarize_set <- function(df, label) {
  z <- df[df$tested_by_edgeR, , drop = FALSE]
  n <- nrow(z)
  k <- sum(z$same_direction, na.rm = TRUE)
  p_binom <- if (n > 0) binom.test(k, n, p = 0.5, alternative = "greater")$p.value else NA_real_
  rho <- if (n >= 3) suppressWarnings(cor(z$primary_logFC, z$logFC, method = "spearman")) else NA_real_
  rho_p <- if (n >= 3) suppressWarnings(cor.test(z$primary_logFC, z$logFC, method = "spearman", exact = FALSE)$p.value) else NA_real_
  data.frame(
    gene_set = label,
    n_prespecified = nrow(df),
    n_represented = sum(df$represented_in_GSE221156),
    n_tested_edgeR = n,
    n_same_direction = k,
    fraction_same_direction = if (n) k/n else NA_real_,
    exact_binomial_one_sided_P = p_binom,
    n_validation_P_lt_0.05 = sum(z$PValue < 0.05, na.rm = TRUE),
    n_validation_FDR_lt_0.05 = sum(z$FDR < 0.05, na.rm = TRUE),
    spearman_rho_original_vs_validation_logFC = rho,
    spearman_P = rho_p,
    stringsAsFactors = FALSE
  )
}

sum27 <- summarize_set(val, "Primary_27_gene_receiver_program")
sum10 <- summarize_set(val[val$is_high_confidence_core, , drop = FALSE], "High_confidence_10_gene_core")
program_summary <- rbind(sum27, sum10)

# Competitive rank-based enrichment of the originally T2D-UP receiver genes and core.
# The two originally down genes are audited individually rather than treated as a 2-gene set.
signed_stat <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
names(signed_stat) <- tab$gene
up27 <- receiver$gene[receiver$primary_logFC > 0]
up27_idx <- which(names(signed_stat) %in% up27)
core_idx <- which(names(signed_stat) %in% core_genes)

camera_rows <- list()
if (length(up27_idx) >= 5) {
  ca <- cameraPR(signed_stat, index = up27_idx, use.ranks = TRUE, inter.gene.cor = 0.01, sort = FALSE)
  camera_rows[[length(camera_rows)+1]] <- data.frame(
    gene_set = "Original_receiver_T2D_up_genes",
    n_prespecified = length(up27), n_tested = length(up27_idx),
    Direction = as.character(ca$Direction[1]), PValue = ca$PValue[1], stringsAsFactors = FALSE)
}
if (length(core_idx) >= 5) {
  ca <- cameraPR(signed_stat, index = core_idx, use.ranks = TRUE, inter.gene.cor = 0.01, sort = FALSE)
  camera_rows[[length(camera_rows)+1]] <- data.frame(
    gene_set = "High_confidence_10_gene_core",
    n_prespecified = length(core_genes), n_tested = length(core_idx),
    Direction = as.character(ca$Direction[1]), PValue = ca$PValue[1], stringsAsFactors = FALSE)
}
camera_res <- if (length(camera_rows)) do.call(rbind, camera_rows) else data.frame()
if (nrow(camera_res)) camera_res$FDR <- p.adjust(camera_res$PValue, method = "BH")

# Individual audit of originally T2D-down genes.
down_audit <- val[val$primary_logFC < 0, , drop = FALSE]

cat("\nINDEPENDENT RECEIVER-PROGRAM VALIDATION\n")
print(program_summary, row.names = FALSE)
cat("\nCompetitive CAMERA-PR enrichment:\n")
if (nrow(camera_res)) print(camera_res, row.names = FALSE) else cat("  Not enough represented genes for CAMERA-PR.\n")
cat("\nOriginally T2D-down receiver genes:\n")
if (nrow(down_audit)) {
  print(down_audit[, c("gene", "primary_logFC", "logFC", "PValue", "FDR", "same_direction")], row.names = FALSE)
} else cat("  None.\n")

write.csv(program_summary,
          file.path(outdir, "Table_P3_GSE221156_ReceiverProgram_ValidationSummary.csv"), row.names = FALSE)
write.csv(camera_res,
          file.path(outdir, "Table_P3_GSE221156_ReceiverProgram_CAMERA_PR.csv"), row.names = FALSE)
write.csv(down_audit,
          file.path(outdir, "Table_P3_GSE221156_ReceiverProgram_DownGeneAudit.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Sensitivity model adjusted for sex + age where estimable
# -----------------------------------------------------------------------------
parse_age <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", x)))
  out[!grepl("^[0-9]+", x)] <- NA_real_
  out
}

sens_result <- NULL
meta_s <- meta_p
meta_s$age_years <- if ("development_stage" %in% names(meta_s)) parse_age(meta_s$development_stage) else NA_real_
meta_s$sex <- if ("sex" %in% names(meta_s)) factor(meta_s$sex) else factor(NA)
complete_s <- complete.cases(meta_s[, c("group", "sex", "age_years")])
meta_s <- meta_s[complete_s, , drop = FALSE]

if (nrow(meta_s) >= 20 && length(unique(meta_s$sex)) >= 2) {
  y_s <- DGEList(counts = counts[, meta_s$donor_id, drop = FALSE])
  keep_s <- filterByExpr(y_s, group = meta_s$group)
  y_s <- y_s[keep_s, , keep.lib.sizes = FALSE]
  y_s <- calcNormFactors(y_s)
  design_s <- model.matrix(~ sex + age_years + group, data = meta_s)
  if (qr(design_s)$rank == ncol(design_s) && nrow(design_s) - ncol(design_s) >= 5) {
    y_s <- estimateDisp(y_s, design_s, robust = TRUE)
    fit_s <- glmQLFit(y_s, design_s, robust = TRUE)
    coef_s <- grep("^groupT2D$", colnames(design_s), value = TRUE)
    if (length(coef_s) == 1) {
      qlf_s <- glmQLFTest(fit_s, coef = coef_s)
      ts <- topTags(qlf_s, n = Inf, sort.by = "none")$table
      ts$gene <- rownames(ts)
      sens_result <- merge(receiver[, c("gene", "primary_logFC")],
                           ts[, c("gene", "logFC", "PValue", "FDR")],
                           by = "gene", all.x = TRUE, sort = FALSE,
                           suffixes = c("_original", "_adjusted"))
      sens_result$same_direction_adjusted <- with(sens_result,
                                                   ifelse(is.na(logFC), NA,
                                                          sign(primary_logFC) == sign(logFC)))
      write.csv(sens_result,
                file.path(outdir, "Table_P3_GSE221156_27Gene_SexAgeAdjusted_Sensitivity.csv"),
                row.names = FALSE)
      cat("\nSEX + AGE ADJUSTED SENSITIVITY\n")
      cat(sprintf("  Donors retained: %d\n", nrow(meta_s)))
      cat(sprintf("  Receiver genes tested: %d\n", sum(!is.na(sens_result$logFC))))
      cat(sprintf("  Same direction: %d/%d\n",
                  sum(sens_result$same_direction_adjusted, na.rm = TRUE),
                  sum(!is.na(sens_result$same_direction_adjusted))))
    }
  } else {
    cat("\nSEX + AGE adjusted sensitivity skipped: design not estimable/full rank.\n")
  }
} else {
  cat("\nSEX + AGE adjusted sensitivity skipped: insufficient complete metadata.\n")
}

# -----------------------------------------------------------------------------
# 6. Independent beta-cell receptor-expression support across donors
# -----------------------------------------------------------------------------
receptor_map <- list(
  RPH3A = c("NRXN1"),
  EDIL3 = c("ITGAV"),
  SPP1 = c("CD44", "ITGAV", "ITGB1", "ITGB5"),
  SEMA7A = c("ITGA1", "ITGB1", "PLXNC1"),
  PPIA = c("BSG"),
  FGF13 = c("FGFR1"),
  FGF12 = c("FGFR1"),
  SST = c("SSTR1", "SSTR2", "SSTR3"),
  CNTN2 = c("CNTN1", "NRCAM", "NRP1"),
  LGALS8 = c("ITGA3"),
  RSPO2 = c("LGR4", "RNF43"),
  FLRT3 = c("UNC5B", "UNC5D"),
  PGF = c("NRP1")
)

# TMM-normalized CPM from the all-48-donor pseudobulk for descriptive support.
y48 <- DGEList(counts = counts[, meta$donor_id, drop = FALSE])
y48 <- calcNormFactors(y48)
cpm48 <- cpm(y48, normalized.lib.sizes = TRUE)

rec_rows <- list()
for (lig in names(receptor_map)) {
  for (rec in receptor_map[[lig]]) {
    if (!rec %in% rownames(cpm48)) {
      rec_rows[[length(rec_rows)+1]] <- data.frame(
        ligand = lig, receptor = rec, present = FALSE,
        disease = NA_character_, n_donors = NA_integer_,
        n_CPM_GE1 = NA_integer_, fraction_CPM_GE1 = NA_real_,
        median_CPM = NA_real_, stringsAsFactors = FALSE)
      next
    }
    for (dis in c("normal", "prediabetes syndrome", "type 2 diabetes mellitus")) {
      ix <- meta$disease == dis
      vv <- cpm48[rec, ix]
      rec_rows[[length(rec_rows)+1]] <- data.frame(
        ligand = lig, receptor = rec, present = TRUE,
        disease = dis, n_donors = sum(ix),
        n_CPM_GE1 = sum(vv >= 1), fraction_CPM_GE1 = mean(vv >= 1),
        median_CPM = median(vv), stringsAsFactors = FALSE)
    }
  }
}
rec_support <- do.call(rbind, rec_rows)
write.csv(rec_support,
          file.path(outdir, "Table_P3_GSE221156_Beta_Receptor_DonorSupport.csv"), row.names = FALSE)

# Compact ligand-level reproduction: receptor supported in >=50% donors in BOTH normal and T2D.
rec_present <- rec_support[rec_support$present & rec_support$disease %in% c("normal", "type 2 diabetes mellitus"), ]
lig_summary <- do.call(rbind, lapply(names(receptor_map), function(lig) {
  x <- rec_present[rec_present$ligand == lig, ]
  receptors <- unique(x$receptor)
  both <- receptors[vapply(receptors, function(r) {
    z <- x[x$receptor == r, ]
    all(c("normal", "type 2 diabetes mellitus") %in% z$disease) && all(z$fraction_CPM_GE1 >= 0.5)
  }, logical(1))]
  data.frame(ligand = lig,
             n_prespecified_receptors = length(receptor_map[[lig]]),
             n_receptors_supported_both = length(both),
             receptors_supported_both = paste(both, collapse = ";"),
             beta_receptor_support_reproduced = length(both) > 0,
             stringsAsFactors = FALSE)
}))
write.csv(lig_summary,
          file.path(outdir, "Table_P3_GSE221156_Ligand_Receptor_Reproduction.csv"), row.names = FALSE)

cat("\nBETA-CELL RECEPTOR SUPPORT IN INDEPENDENT 48-DONOR COHORT\n")
cat(sprintf("  Candidate ligands with >=1 receptor at CPM>=1 in >=50%% of BOTH normal and T2D donors: %d/%d\n",
            sum(lig_summary$beta_receptor_support_reproduced), nrow(lig_summary)))
print(lig_summary, row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Secondary descriptive normal -> prediabetes -> T2D receiver score
# -----------------------------------------------------------------------------
# Direction-align each receiver gene to the original GSE81608 effect, z-score across donors,
# then average per donor. This is descriptive/secondary, not the primary validation test.
receiver_present <- receiver$gene[receiver$gene %in% rownames(cpm48)]
logcpm48 <- log2(cpm48 + 0.25)
if (length(receiver_present) >= 5) {
  mat <- logcpm48[receiver_present, , drop = FALSE]
  zmat <- t(scale(t(mat)))
  orig_sign <- sign(receiver$primary_logFC[match(receiver_present, receiver$gene)])
  zmat <- zmat * orig_sign
  score <- colMeans(zmat, na.rm = TRUE)
  trajectory <- data.frame(
    donor_id = colnames(zmat),
    disease = meta$disease[match(colnames(zmat), meta$donor_id)],
    receiver_direction_aligned_score = score,
    stringsAsFactors = FALSE
  )
  write.csv(trajectory,
            file.path(outdir, "Table_P3_GSE221156_ReceiverScore_All48Donors.csv"), row.names = FALSE)
  cat("\nSECONDARY DIRECTION-ALIGNED RECEIVER SCORE (descriptive)\n")
  med <- aggregate(receiver_direction_aligned_score ~ disease, trajectory, median)
  print(med, row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 8. Save reproducibility object and session information
# -----------------------------------------------------------------------------
saveRDS(list(
  primary_edgeR = tab,
  receiver_validation = val,
  receiver_program_summary = program_summary,
  cameraPR = camera_res,
  receptor_support = rec_support,
  ligand_receptor_reproduction = lig_summary,
  primary_metadata = meta_p,
  filter_keep = keep_gene
), file.path("data/processed", "GSE221156_Beta_IndependentValidation_13d.rds"))

capture.output(sessionInfo(), file = file.path(reprodir, "13d_GSE221156_validation_sessionInfo.txt"))

cat("\nOUTPUT FILES\n")
for (f in c(
  de_file,
  receiver_out,
  file.path(outdir, "Table_P3_GSE221156_ReceiverProgram_ValidationSummary.csv"),
  file.path(outdir, "Table_P3_GSE221156_ReceiverProgram_CAMERA_PR.csv"),
  file.path(outdir, "Table_P3_GSE221156_ReceiverProgram_DownGeneAudit.csv"),
  file.path(outdir, "Table_P3_GSE221156_Beta_Receptor_DonorSupport.csv"),
  file.path(outdir, "Table_P3_GSE221156_Ligand_Receptor_Reproduction.csv"),
  file.path(outdir, "Table_P3_GSE221156_ReceiverScore_All48Donors.csv"),
  file.path("data/processed", "GSE221156_Beta_IndependentValidation_13d.rds"),
  file.path(reprodir, "13d_GSE221156_validation_sessionInfo.txt")
)) if (file.exists(f)) cat("  ", f, "\n", sep = "")

cat("\nINTERPRETATION GUARDRAIL\n")
cat("  This is independent beta-cell transcriptional/receptor validation.\n")
cat("  It does NOT demonstrate direct sensory-neuron-to-beta-cell signaling or causality.\n")
cat("  Prediabetes trajectory scores are secondary/descriptive.\n")

cat("\n", strrep("=", 78), "\n", sep = "")
cat("SCRIPT 13d COMPLETE\n")
cat(strrep("=", 78), "\n", sep = "")
