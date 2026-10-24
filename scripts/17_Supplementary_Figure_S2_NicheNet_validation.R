#!/usr/bin/env Rscript

# Supplementary Figure S2 — Independent NicheNet replication and empirical-null sensitivity
# Paper 3: Integrated human transcriptomics prioritizes candidate sensory neuron–islet signaling axes in T2D
#
# This script reads only frozen outputs from Scripts 16a, 16b and 16c.
# It does NOT rerun NicheNet, differential expression, or permutations.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

cat(strrep("=", 78), "\n", sep = "")
cat("SUPPLEMENTARY FIGURE S2 — NICHE-NET VALIDATION / EMPIRICAL NULL\n")
cat(strrep("=", 78), "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
root <- getwd()

step1_file <- file.path(
  root,
  "results/validation/GSE221156/NicheNet_replication",
  "Table_P3_GSE221156_vs_GSE81608_NicheNet_RankReplication.csv"
)

step2_file <- file.path(
  root,
  "results/validation/GSE221156/NicheNet_broad_replication",
  "Table_P3_BroadHC_NicheNet_GSE81608_vs_GSE221156.csv"
)

perm_dir <- file.path(
  root,
  "results/validation/GSE221156/NicheNet_permutation/final"
)

out_dir <- file.path(root, "figures/supplementary")
src_dir <- file.path(root, "results/figure_source_data/Supplementary_Figure_S2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(step1_file, step2_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing frozen input file(s):\n", paste(" -", missing_files, collapse = "\n"))
}
if (!dir.exists(perm_dir)) stop("Permutation final directory not found: ", perm_dir)

cat("Frozen Step 1 file: PASS\n")
cat("Frozen Step 2 file: PASS\n")
cat("Frozen Step 3 directory: PASS\n")

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
read_csv_safe <- function(f) {
  tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
}

has_cols <- function(x, cols) {
  !is.null(x) && all(cols %in% names(x))
}

main_genes <- c("SEMA7A", "EDIL3", "FLRT3", "PGF", "SPP1")
main_order <- c("SEMA7A", "EDIL3", "FLRT3", "PGF", "SPP1")

# Publication theme: restrained, readable in grayscale as well as color.
theme_nc <- theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
    plot.subtitle = element_text(size = 9.2, hjust = 0, margin = margin(b = 6)),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(8, 10, 8, 8)
  )

# -----------------------------------------------------------------------------
# 3. Step 1 — 10-candidate discovery -> independent receiver rank replication
# -----------------------------------------------------------------------------
rank10 <- read.csv(step1_file, check.names = FALSE)
need1 <- c("test_ligand", "discovery_rank", "validation_rank")
if (!all(need1 %in% names(rank10))) {
  stop("Step 1 rank-replication file lacks required columns: ", paste(need1, collapse = ", "))
}
if (nrow(rank10) != 10) stop("Expected exactly 10 NicheNet-evaluable sensory candidates in Step 1; found ", nrow(rank10))

rank10$test_ligand <- as.character(rank10$test_ligand)
rank10 <- rank10[order(rank10$discovery_rank), , drop = FALSE]
rank10$ligand_factor <- factor(rank10$test_ligand, levels = rev(rank10$test_ligand))

rank10_long <- rbind(
  data.frame(test_ligand = rank10$test_ligand, ligand_factor = rank10$ligand_factor,
             cohort = "GSE81608 discovery", rank = rank10$discovery_rank),
  data.frame(test_ligand = rank10$test_ligand, ligand_factor = rank10$ligand_factor,
             cohort = "GSE221156 independent", rank = rank10$validation_rank)
)
rank10_long$cohort <- factor(rank10_long$cohort,
                             levels = c("GSE81608 discovery", "GSE221156 independent"))

pA <- ggplot(rank10, aes(y = ligand_factor)) +
  geom_segment(aes(x = discovery_rank, xend = validation_rank,
                   yend = ligand_factor), linewidth = 0.55, color = "grey65") +
  geom_point(aes(x = discovery_rank), shape = 21, size = 3.0,
             stroke = 0.65, fill = "white") +
  geom_point(aes(x = validation_rank), shape = 16, size = 2.7) +
  geom_text(aes(x = discovery_rank, label = discovery_rank),
            nudge_y = 0.20, size = 2.6, vjust = 0) +
  geom_text(aes(x = validation_rank, label = validation_rank),
            nudge_y = -0.20, size = 2.6, vjust = 1) +
  scale_x_continuous(breaks = 1:10, limits = c(0.5, 10.5)) +
  labs(
    title = "A  Sensory-candidate rank replication",
    subtitle = "Open = GSE81608; filled = independent GSE221156 receiver",
    x = "NicheNet corrected-AUPR rank (1 = highest)", y = NULL
  ) +
  theme_nc

# -----------------------------------------------------------------------------
# 4. Step 2 — broad 42-ligand cross-cohort rank concordance
# -----------------------------------------------------------------------------
broad <- read.csv(step2_file, check.names = FALSE)
need2 <- c("test_ligand", "discovery_broad_rank", "validation_broad_rank")
if (!all(need2 %in% names(broad))) {
  stop("Step 2 broad-replication file lacks required columns: ", paste(need2, collapse = ", "))
}
if (nrow(broad) != 42) stop("Expected exactly 42 NicheNet-evaluable broad-universe ligands; found ", nrow(broad))

broad$test_ligand <- as.character(broad$test_ligand)
broad$is_main <- broad$test_ligand %in% main_genes
rho_broad <- suppressWarnings(cor(broad$discovery_broad_rank, broad$validation_broad_rank,
                                  method = "spearman", use = "complete.obs"))
p_broad <- suppressWarnings(cor.test(broad$discovery_broad_rank, broad$validation_broad_rank,
                                     method = "spearman", exact = FALSE)$p.value)

main_broad <- broad[broad$is_main, , drop = FALSE]
# Small deterministic label nudges to avoid overlap without an extra package.
main_broad$nudge_x <- c(SEMA7A = -2.0, EDIL3 = -2.0, FLRT3 = 2.2, PGF = 2.0, SPP1 = -2.2)[main_broad$test_ligand]
main_broad$nudge_y <- c(SEMA7A = -2.0, EDIL3 = 2.0, FLRT3 = -2.2, PGF = 2.0, SPP1 = -2.0)[main_broad$test_ligand]

pB <- ggplot(broad, aes(discovery_broad_rank, validation_broad_rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.45, color = "grey55") +
  geom_point(aes(shape = is_main), size = 2.5, stroke = 0.65) +
  geom_text(
    data = main_broad,
    aes(x = discovery_broad_rank + nudge_x,
        y = validation_broad_rank + nudge_y,
        label = test_ligand),
    fontface = "bold", size = 3.0, check_overlap = FALSE
  ) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  scale_x_reverse(breaks = c(40, 30, 20, 10, 1), limits = c(42.5, 0.5)) +
  scale_y_reverse(breaks = c(40, 30, 20, 10, 1), limits = c(42.5, 0.5)) +
  coord_fixed() +
  annotate("text", x = 39, y = 3.0,
           label = sprintf("Spearman rho = %.3f\nP = %.2e", rho_broad, p_broad),
           hjust = 0, vjust = 0, size = 3.0) +
  labs(
    title = "B  Broad-ligand rank concordance",
    subtitle = "Prespecified 42-ligand high-confidence beta-supported universe",
    x = "GSE81608 broad rank", y = "GSE221156 broad rank"
  ) +
  theme_nc

# -----------------------------------------------------------------------------
# 5. Step 3 — discover the final 10,000-permutation summary table(s)
# -----------------------------------------------------------------------------
perm_csvs <- list.files(perm_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
if (length(perm_csvs) == 0) stop("No CSV outputs found in final permutation directory: ", perm_dir)
perm_tabs <- setNames(lapply(perm_csvs, read_csv_safe), basename(perm_csvs))

# Preferred: one cross-cohort table already containing both cohorts.
combined_required <- c(
  "test_ligand",
  "empirical_P_GSE81608", "empirical_FDR_GSE81608",
  "empirical_P_GSE221156", "empirical_FDR_GSE221156"
)
combined_idx <- which(vapply(perm_tabs, has_cols, logical(1), cols = combined_required))

if (length(combined_idx) >= 1) {
  perm <- perm_tabs[[combined_idx[1]]]
  perm_source <- names(perm_tabs)[combined_idx[1]]
} else {
  # Fallback: merge cohort-specific empirical-null tables.
  single_required <- c("test_ligand", "empirical_P", "empirical_FDR")
  single_idx <- which(vapply(perm_tabs, has_cols, logical(1), cols = single_required))
  if (length(single_idx) < 2) {
    stop(
      "Could not identify final permutation summary table(s).\n",
      "Expected either a combined table with columns:\n  ",
      paste(combined_required, collapse = ", "),
      "\nor two cohort-specific tables with: test_ligand, empirical_P, empirical_FDR.\n",
      "CSV files found:\n  ", paste(basename(perm_csvs), collapse = "\n  ")
    )
  }

  candidates <- perm_tabs[single_idx]
  get_cohort <- function(tab, nm) {
    if ("cohort" %in% names(tab) && length(unique(tab$cohort)) == 1) return(as.character(unique(tab$cohort)))
    if (grepl("81608", nm, ignore.case = TRUE)) return("GSE81608")
    if (grepl("221156", nm, ignore.case = TRUE)) return("GSE221156")
    NA_character_
  }
  cohort_names <- mapply(get_cohort, candidates, names(candidates), SIMPLIFY = TRUE)
  i816 <- which(cohort_names == "GSE81608")[1]
  i221 <- which(cohort_names == "GSE221156")[1]
  if (is.na(i816) || is.na(i221)) stop("Found cohort-specific tables but could not identify both GSE81608 and GSE221156.")

  a <- candidates[[i816]][, c("test_ligand", "empirical_P", "empirical_FDR")]
  b <- candidates[[i221]][, c("test_ligand", "empirical_P", "empirical_FDR")]
  names(a)[2:3] <- c("empirical_P_GSE81608", "empirical_FDR_GSE81608")
  names(b)[2:3] <- c("empirical_P_GSE221156", "empirical_FDR_GSE221156")
  perm <- merge(a, b, by = "test_ligand", all = FALSE, sort = FALSE)
  perm_source <- paste(names(candidates)[c(i816, i221)], collapse = " + ")
}

perm$test_ligand <- as.character(perm$test_ligand)
if (nrow(perm) != 42) stop("Expected 42 ligands in final permutation summary; found ", nrow(perm))

# Hard audit for final run, preventing accidental use of the 20-permutation smoke test.
min_possible <- 1 / 10001
if (min(perm$empirical_P_GSE221156, na.rm = TRUE) > min_possible * 1.01) {
  stop("Permutation table does not appear to be the final 10,000-permutation run. Minimum empirical P = ",
       signif(min(perm$empirical_P_GSE221156, na.rm = TRUE), 6),
       "; expected minimum near 1/10001 = ", signif(min_possible, 6))
}

n_sig_816 <- sum(perm$empirical_FDR_GSE81608 < 0.05, na.rm = TRUE)
n_sig_221 <- sum(perm$empirical_FDR_GSE221156 < 0.05, na.rm = TRUE)
n_sig_both <- sum(perm$empirical_FDR_GSE81608 < 0.05 & perm$empirical_FDR_GSE221156 < 0.05, na.rm = TRUE)

# Frozen final expectations from Script 16c.
if (n_sig_816 != 0) stop("Audit failed: expected 0/42 GSE81608 ligands at empirical FDR < 0.05; found ", n_sig_816)
if (n_sig_221 != 41) stop("Audit failed: expected 41/42 GSE221156 ligands at empirical FDR < 0.05; found ", n_sig_221)
if (n_sig_both != 0) stop("Audit failed: expected 0/42 ligands significant in both cohorts; found ", n_sig_both)

cat("Final permutation table: ", perm_source, "\n", sep = "")
cat("Final permutation audit: PASS (0/42 GSE81608; 41/42 GSE221156; 0/42 both)\n")

# -----------------------------------------------------------------------------
# 6. Panel C — five principal candidates, empirical FDR by cohort
# -----------------------------------------------------------------------------
main_perm <- perm[match(main_order, perm$test_ligand), , drop = FALSE]
if (any(is.na(main_perm$test_ligand))) stop("One or more main candidates missing from permutation table.")

main_long <- rbind(
  data.frame(test_ligand = main_perm$test_ligand,
             cohort = "GSE81608 discovery",
             empirical_FDR = main_perm$empirical_FDR_GSE81608,
             empirical_P = main_perm$empirical_P_GSE81608),
  data.frame(test_ligand = main_perm$test_ligand,
             cohort = "GSE221156 independent",
             empirical_FDR = main_perm$empirical_FDR_GSE221156,
             empirical_P = main_perm$empirical_P_GSE221156)
)
main_long$test_ligand <- factor(main_long$test_ligand, levels = rev(main_order))
main_long$cohort <- factor(main_long$cohort,
                           levels = c("GSE81608 discovery", "GSE221156 independent"))
main_long$minuslog10FDR <- -log10(pmax(main_long$empirical_FDR, .Machine$double.xmin))

pC <- ggplot(main_long, aes(x = minuslog10FDR, y = test_ligand, shape = cohort)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.5, color = "grey45") +
  geom_point(size = 3.1, stroke = 0.75) +
  scale_shape_manual(values = c("GSE81608 discovery" = 1, "GSE221156 independent" = 16)) +
  annotate("text", x = -log10(0.05), y = 5.45, label = "FDR = 0.05",
           hjust = -0.05, vjust = 1, size = 2.8) +
  labs(
    title = "C  Expression-matched empirical null",
    subtitle = "Five principal candidates; 10,000 permutations per cohort",
    x = expression(-log[10]("empirical FDR")), y = NULL
  ) +
  theme_nc +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# 7. Panel D — all 42 ligands, empirical-FDR landscape
# -----------------------------------------------------------------------------
perm$is_main <- perm$test_ligand %in% main_genes
perm$xFDR <- -log10(pmax(perm$empirical_FDR_GSE81608, .Machine$double.xmin))
perm$yFDR <- -log10(pmax(perm$empirical_FDR_GSE221156, .Machine$double.xmin))
main_fdr <- perm[perm$is_main, , drop = FALSE]

# Deterministic offsets for five labels.
main_fdr$dx <- c(SEMA7A = 0.05, EDIL3 = 0.05, FLRT3 = 0.05, PGF = 0.05, SPP1 = 0.05)[main_fdr$test_ligand]
main_fdr$dy <- c(SEMA7A = 0.14, EDIL3 = -0.18, FLRT3 = 0.28, PGF = -0.32, SPP1 = 0.00)[main_fdr$test_ligand]

thr <- -log10(0.05)
pD <- ggplot(perm, aes(xFDR, yFDR)) +
  geom_vline(xintercept = thr, linetype = 2, linewidth = 0.5, color = "grey50") +
  geom_hline(yintercept = thr, linetype = 2, linewidth = 0.5, color = "grey50") +
  geom_point(aes(shape = is_main), size = 2.5, stroke = 0.65) +
  geom_text(
    data = main_fdr,
    aes(x = xFDR + dx, y = yFDR + dy, label = test_ligand),
    fontface = "bold", hjust = 0, size = 2.9
  ) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("FDR < 0.05\nGSE81608: %d/42\nGSE221156: %d/42\nBoth: %d/42",
                           n_sig_816, n_sig_221, n_sig_both),
           hjust = 1.05, vjust = -0.10, size = 2.8, label.size = 0.25) +
  labs(
    title = "D  Empirical-null specificity across 42 ligands",
    subtitle = "No ligand reaches empirical FDR < 0.05 in both cohorts",
    x = expression(GSE81608~~-log[10]("empirical FDR")),
    y = expression(GSE221156~~-log[10]("empirical FDR"))
  ) +
  theme_nc

# -----------------------------------------------------------------------------
# 8. Assemble and save
# -----------------------------------------------------------------------------
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Independent NicheNet replication and expression-matched empirical-null sensitivity",
    subtitle = paste0(
      "Rank-based replication is shown separately from empirical significance; ",
      "NicheNet activity indicates regulatory compatibility, not direct sensory neuron–beta-cell signaling."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 9.5, hjust = 0)
    )
  )

# Source-data exports used by the figure.
write.csv(rank10, file.path(src_dir, "Figure_S2_PanelA_SensoryCandidateRankReplication.csv"), row.names = FALSE)
write.csv(broad, file.path(src_dir, "Figure_S2_PanelB_BroadLigandRankConcordance.csv"), row.names = FALSE)
write.csv(main_long, file.path(src_dir, "Figure_S2_PanelC_MainCandidateEmpiricalFDR.csv"), row.names = FALSE)
write.csv(perm, file.path(src_dir, "Figure_S2_PanelD_AllLigandEmpiricalFDR.csv"), row.names = FALSE)

pdf_file <- file.path(out_dir, "Supplementary_Figure_S2_NicheNet_Validation_EmpiricalNull.pdf")
png_file <- file.path(out_dir, "Supplementary_Figure_S2_NicheNet_Validation_EmpiricalNull.png")
tif_file <- file.path(out_dir, "Supplementary_Figure_S2_NicheNet_Validation_EmpiricalNull.tiff")

ggsave(pdf_file, combined, width = 12.6, height = 9.4, units = "in", device = cairo_pdf)
ggsave(png_file, combined, width = 12.6, height = 9.4, units = "in", dpi = 450)
ggsave(tif_file, combined, width = 12.6, height = 9.4, units = "in", dpi = 600, compression = "lzw")

# Save individual panels too.
for (nm in names(list(A = pA, B = pB, C = pC, D = pD))) {
  pp <- list(A = pA, B = pB, C = pC, D = pD)[[nm]]
  ggsave(file.path(out_dir, paste0("Supplementary_Figure_S2_Panel", nm, ".pdf")),
         pp, width = 6.2, height = 4.5, units = "in", device = cairo_pdf)
}

# Reproducibility record.
record <- c(
  paste0("Generated: ", format(Sys.time())),
  paste0("Step 1 input: ", step1_file),
  paste0("Step 2 input: ", step2_file),
  paste0("Step 3 source: ", perm_source),
  paste0("Broad rank Spearman rho: ", signif(rho_broad, 7)),
  paste0("Broad rank Spearman P: ", signif(p_broad, 7)),
  paste0("Empirical FDR < 0.05 GSE81608: ", n_sig_816, "/42"),
  paste0("Empirical FDR < 0.05 GSE221156: ", n_sig_221, "/42"),
  paste0("Empirical FDR < 0.05 both: ", n_sig_both, "/42"),
  "",
  capture.output(sessionInfo())
)
writeLines(record, file.path(src_dir, "Supplementary_Figure_S2_sessionInfo.txt"))

cat("\nAUDIT\n")
cat("Step 1 sensory candidates: ", nrow(rank10), "\n", sep = "")
cat("Step 2 broad ligands: ", nrow(broad), "\n", sep = "")
cat(sprintf("Broad rank Spearman rho = %.5f; P = %.5g\n", rho_broad, p_broad))
cat("Step 3 GSE81608 empirical FDR < 0.05: ", n_sig_816, "/42\n", sep = "")
cat("Step 3 GSE221156 empirical FDR < 0.05: ", n_sig_221, "/42\n", sep = "")
cat("Step 3 BOTH empirical FDR < 0.05: ", n_sig_both, "/42\n", sep = "")
cat("\nOUTPUTS\n")
cat(" - ", pdf_file, "\n", sep = "")
cat(" - ", png_file, "\n", sep = "")
cat(" - ", tif_file, "\n", sep = "")
cat(" - source data: ", src_dir, "\n", sep = "")
cat("\nINTERPRETATION GUARDRAIL\n")
cat(" Rank replication and empirical significance are distinct evidence layers.\n")
cat(" The 42-ligand universe and candidate set were fixed before permutation testing.\n")
cat(" Empirical-null significance assesses receiver-program specificity within the frozen NicheNet prior.\n")
cat(" It does not demonstrate direct inter-organ signaling or causality.\n")
cat("\n", strrep("=", 78), "\n", sep = "")
cat("SUPPLEMENTARY FIGURE S2 COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 78), "\n", sep = "")
