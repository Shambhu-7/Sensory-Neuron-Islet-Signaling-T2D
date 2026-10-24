#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 10
# DRG ligand discovery and human islet receptor compatibility
# Figure-only script using frozen results from Scripts 03/04.
# No ligand discovery, receptor analysis, or statistics are rerun.
# ============================================================

options(stringsAsFactors = FALSE)

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

cat("====================================================\n")
cat("SCRIPT 10: DRG LIGAND DISCOVERY + ISLET COMPATIBILITY FIGURE\n")
cat("Frozen-result visualization only; no analysis is rerun.\n")
cat("====================================================\n\n")

# -----------------------------
# Paths
# -----------------------------
f_all_ligands <- "results/ligands/Table_P3_All_Curated_Ligands_with_DE_Robustness.csv"
f_prioritized <- "results/ligands/Table_P3_Prioritized_Robust_DRG_Ligands.csv"
f_shortlist <- "results/islet/Table_P3_SensorySource_BetaCell_HC_Shortlist.csv"
f_lr_sens <- "results/islet/Table_P3_HC_LR_1492_vs_1600_Sensitivity.csv"
f_lig_sens <- "results/islet/Table_P3_HC_Ligand_1492_vs_1600_Sensitivity.csv"

required_files <- c(f_all_ligands, f_prioritized, f_shortlist, f_lr_sens, f_lig_sens)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required frozen input file(s):\n", paste(missing_files, collapse = "\n"))
}

dir.create("figures/ligands", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figure_source_data", recursive = TRUE, showWarnings = FALSE)
dir.create("reproducibility", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read frozen results
# -----------------------------
all_lig <- read.csv(f_all_ligands, check.names = FALSE)
prior <- read.csv(f_prioritized, check.names = FALSE)
short <- read.csv(f_shortlist, check.names = FALSE)
lr_sens <- read.csv(f_lr_sens, check.names = FALSE)
lig_sens <- read.csv(f_lig_sens, check.names = FALSE)

stopifnot(nrow(short) == 13)

# -----------------------------
# Helpers
# -----------------------------
tier_label <- function(x) {
  dplyr::recode(
    x,
    "Tier_A_strict" = "Tier A (strict)",
    "Tier_B_LODO_robust" = "Tier B (LODO-robust)",
    "Tier_C_directionally_robust" = "Tier C (direction-robust)",
    .default = x
  )
}

boolify <- function(x) {
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")
}

base_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9.2, hjust = 0, margin = margin(b = 7)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 8.8),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(7, 8, 7, 8)
  )

# ============================================================
# PANEL A — ligand discovery funnel
# ============================================================

n_tested <- nrow(all_lig)
n_prior <- nrow(prior)

tier_counts <- prior %>%
  count(candidate_tier, name = "n") %>%
  mutate(label = tier_label(candidate_tier))

# Validate the expected frozen analysis totals.
if (n_tested != 283) warning("Expected 283 DE-tested curated ligands; observed ", n_tested)
if (n_prior != 149) warning("Expected 149 prioritized robust ligands; observed ", n_prior)

n_A <- sum(prior$candidate_tier == "Tier_A_strict", na.rm = TRUE)
n_B <- sum(prior$candidate_tier == "Tier_B_LODO_robust", na.rm = TRUE)
n_C <- sum(prior$candidate_tier == "Tier_C_directionally_robust", na.rm = TRUE)

funnel <- data.frame(
  stage = factor(
    c("DE-tested curated LR ligands",
      "Prioritized robust DRG ligands",
      "Tier A (strict)",
      "Tier B (LODO-robust)",
      "Tier C (direction-robust)",
      "Sensory-source + beta-compatible"),
    levels = rev(c("DE-tested curated LR ligands",
                   "Prioritized robust DRG ligands",
                   "Tier A (strict)",
                   "Tier B (LODO-robust)",
                   "Tier C (direction-robust)",
                   "Sensory-source + beta-compatible"))
  ),
  n = c(n_tested, n_prior, n_A, n_B, n_C, nrow(short)),
  group = c("Discovery", "Discovery", "Tier A", "Tier B", "Tier C", "Final shortlist")
)

pA <- ggplot(funnel, aes(x = n, y = stage, fill = group)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.4, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "A  Ligand discovery and prioritization",
    subtitle = "Frozen OmniPath-derived ligand-receptor network and donor-robust DRG filtering",
    x = "Number of ligands",
    y = NULL
  ) +
  base_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ============================================================
# PANEL B — 13 sensory-source candidates
# ============================================================

short_plot <- short %>%
  mutate(
    nearby_enrichment = -logFC,
    neglogFDR = -log10(pmax(FDR, .Machine$double.xmin)),
    tier = tier_label(candidate_tier),
    ligand = factor(ligand, levels = rev(ligand))
  )

pB <- ggplot(short_plot, aes(y = ligand, x = nearby_enrichment)) +
  geom_segment(aes(x = 0, xend = nearby_enrichment, yend = ligand),
               linewidth = 0.55, color = "grey45") +
  geom_point(aes(shape = tier, size = neglogFDR), fill = "black", color = "black") +
  scale_shape_manual(values = c(
    "Tier A (strict)" = 16,
    "Tier B (LODO-robust)" = 17,
    "Tier C (direction-robust)" = 15
  )) +
  scale_size_continuous(range = c(2.2, 5.2), breaks = c(2, 3, 4, 5, 6)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.08))) +
  labs(
    title = "B  High-confidence LR-filtered sensory-source candidates",
    subtitle = "All 13 are nearby-neuron-associated with β-cell receptor compatibility",
    x = expression("Nearby-neuron enrichment (" * -log[2] * "FC)"),
    y = NULL,
    shape = "DRG tier",
    size = expression(-log[10] * " FDR")
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 0)
  ) +
  guides(
    shape = guide_legend(order = 1, nrow = 1),
    size = guide_legend(order = 2, nrow = 1)
  )

# ============================================================
# PANEL C — beta-cell receptor compatibility for 13 candidates
# ============================================================

receptor_plot <- short %>%
  transmute(
    ligand,
    Supported = n_HC_beta_both_receptors,
    Strong = n_HC_beta_strong_both_receptors
  ) %>%
  pivot_longer(c(Supported, Strong), names_to = "support", values_to = "n_receptors") %>%
  mutate(
    support = factor(
      support,
      levels = c("Supported", "Strong"),
      labels = c("Supported in ND + T2D", "Strong in ND + T2D")
    ),
    ligand = factor(ligand, levels = rev(short$ligand))
  )

pC <- ggplot(receptor_plot, aes(x = n_receptors, y = ligand, fill = support)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(
    aes(label = ifelse(n_receptors > 0, n_receptors, "")),
    position = position_dodge(width = 0.72),
    hjust = -0.2,
    size = 2.8
  ) +
  scale_x_continuous(breaks = 0:4, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "C  Human β-cell receptor compatibility",
    subtitle = "High-confidence receptors supported in both non-diabetic and T2D donors",
    x = "Number of β-cell receptors",
    y = NULL,
    fill = NULL
  ) +
  base_theme +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

# ============================================================
# PANEL D — independent 1600-cell sensitivity representation
# ============================================================

# Broader high-confidence beta-supported pool:
# derive the summary directly from frozen sensitivity tables.
lig_sens2 <- lig_sens %>%
  mutate(
    primary1492_beta_supported = boolify(primary1492_beta_supported),
    sens1600_beta_supported = boolify(sens1600_beta_supported),
    beta_supported_in_both_representations = boolify(beta_supported_in_both_representations)
  )

lr_sens2 <- lr_sens %>%
  mutate(
    primary1492_beta_both = boolify(primary1492_beta_both),
    sens1600_beta_both = boolify(sens1600_beta_both),
    beta_primary_class_agreement = boolify(beta_primary_class_agreement),
    beta_primary_reproduced = boolify(beta_primary_reproduced)
  )

primary_pairs <- sum(lr_sens2$primary1492_beta_both, na.rm = TRUE)
reproduced_pairs <- sum(lr_sens2$primary1492_beta_both & lr_sens2$sens1600_beta_both, na.rm = TRUE)

primary_ligands <- sum(lig_sens2$primary1492_beta_supported, na.rm = TRUE)
reproduced_ligands <- sum(
  lig_sens2$primary1492_beta_supported & lig_sens2$sens1600_beta_supported,
  na.rm = TRUE
)

pair_agreement <- mean(lr_sens2$beta_primary_class_agreement, na.rm = TRUE)

# Shortlist-specific support reproducibility.
near_sens_file <- "results/islet/Table_P3_NearbyNeuron_HC_1492_vs_1600_Sensitivity.csv"
if (!file.exists(near_sens_file)) stop("Missing: ", near_sens_file)
near_sens <- read.csv(near_sens_file, check.names = FALSE) %>%
  mutate(
    beta_supported_in_both_representations =
      boolify(beta_supported_in_both_representations),
    strong_beta_supported_in_both_representations =
      boolify(strong_beta_supported_in_both_representations)
  )

n_short_both <- sum(near_sens$beta_supported_in_both_representations, na.rm = TRUE)
n_short_strong_both <- sum(near_sens$strong_beta_supported_in_both_representations, na.rm = TRUE)

sens_summary <- data.frame(
  metric = factor(
    c("Primary LR pairs reproduced",
      "Primary ligands retained",
      "Pair-level classification agreement",
      "13-candidate β-cell support reproduced",
      "13-candidate strong support reproduced"),
    levels = rev(c("Primary LR pairs reproduced",
                   "Primary ligands retained",
                   "Pair-level classification agreement",
                   "13-candidate β-cell support reproduced",
                   "13-candidate strong support reproduced"))
  ),
  value = c(
    ifelse(primary_pairs > 0, reproduced_pairs / primary_pairs, NA_real_),
    ifelse(primary_ligands > 0, reproduced_ligands / primary_ligands, NA_real_),
    pair_agreement,
    n_short_both / nrow(short),
    n_short_strong_both / nrow(short)
  ),
  label = c(
    paste0(reproduced_pairs, "/", primary_pairs),
    paste0(reproduced_ligands, "/", primary_ligands),
    sprintf("%.1f%%", 100 * pair_agreement),
    paste0(n_short_both, "/", nrow(short)),
    paste0(n_short_strong_both, "/", nrow(short))
  )
)

pD <- ggplot(sens_summary, aes(x = value, y = metric)) +
  geom_segment(aes(x = 0, xend = value, yend = metric),
               linewidth = 1.0, color = "grey70") +
  geom_point(size = 4.2) +
  geom_text(aes(label = label), hjust = -0.35, size = 3.3, fontface = "bold") +
  scale_x_continuous(
    limits = c(0, 1.16),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "D  GSE81608 representation sensitivity",
    subtitle = "Primary 1,492-cell receptor classifications compared with\nthe 1,600-cell GEO representation",
    x = "Reproducibility / agreement",
    y = NULL
  ) +
  base_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ============================================================
# Source-data freeze for the figure
# ============================================================

write.csv(funnel,
          "results/figure_source_data/Figure_P3_LigandIslet_PanelA_counts.csv",
          row.names = FALSE)
write.csv(short_plot,
          "results/figure_source_data/Figure_P3_LigandIslet_PanelB_13candidates.csv",
          row.names = FALSE)
write.csv(receptor_plot,
          "results/figure_source_data/Figure_P3_LigandIslet_PanelC_receptors.csv",
          row.names = FALSE)
write.csv(sens_summary,
          "results/figure_source_data/Figure_P3_LigandIslet_PanelD_sensitivity.csv",
          row.names = FALSE)

# ============================================================
# Combined figure
# ============================================================

combined <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Robust DRG ligand discovery identifies β-cell-compatible sensory-source candidates",
    subtitle = "Donor-robust human diabetic DRG signals integrated with human β-cell receptor compatibility and GSE81608 sensitivity analysis",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 10.5, hjust = 0, margin = margin(b = 8))
    )
  ) +
  plot_layout(heights = c(1, 1.05))

# -----------------------------
# Save individual panels
# -----------------------------
panels <- list(
  A_Discovery = pA,
  B_13Candidates = pB,
  C_BetaReceptors = pC,
  D_Sensitivity = pD
)

for (nm in names(panels)) {
  ggsave(
    filename = paste0("figures/ligands/Figure_P3_LigandIslet_10b_", nm, ".pdf"),
    plot = panels[[nm]],
    device = cairo_pdf,
    width = 7.0, height = 5.1, units = "in"
  )
  ggsave(
    filename = paste0("figures/ligands/Figure_P3_LigandIslet_10b_", nm, ".png"),
    plot = panels[[nm]],
    width = 7.0, height = 5.1, units = "in", dpi = 300
  )
}

# -----------------------------
# Save combined figure
# -----------------------------
ggsave(
  "figures/ligands/Figure_P3_DRG_Ligand_IsletCompatibility_10b.pdf",
  combined,
  device = cairo_pdf,
  width = 14.0, height = 10.0, units = "in"
)

ggsave(
  "figures/ligands/Figure_P3_DRG_Ligand_IsletCompatibility_10b.png",
  combined,
  width = 14.0, height = 10.0, units = "in", dpi = 300
)

ggsave(
  "figures/ligands/Figure_P3_DRG_Ligand_IsletCompatibility_10b.tiff",
  combined,
  width = 14.0, height = 10.0, units = "in",
  dpi = 600, compression = "lzw"
)

# ============================================================
# Audit
# ============================================================

cat("\n====================================================\n")
cat("SCRIPT 10b FIGURE AUDIT\n")
cat("====================================================\n")
cat("DE-tested curated LR ligands:", n_tested, "\n")
cat("Prioritized robust DRG ligands:", n_prior, "\n")
cat("Tier A:", n_A, "\n")
cat("Tier B:", n_B, "\n")
cat("Tier C:", n_C, "\n")
cat("Sensory-source beta-compatible shortlist:", nrow(short), "\n")
cat("Shortlist Tier A:", sum(short$candidate_tier == "Tier_A_strict"), "\n")
cat("Shortlist Tier B:", sum(short$candidate_tier == "Tier_B_LODO_robust"), "\n")
cat("Shortlist Tier C:", sum(short$candidate_tier == "Tier_C_directionally_robust"), "\n")
cat("Shortlist beta-supported:", sum(boolify(short$beta_HC_supported)), "/", nrow(short), "\n")
cat("Shortlist strong beta-supported:", sum(boolify(short$beta_HC_strong)), "/", nrow(short), "\n")
cat("\nSensitivity representation:\n")
cat("Primary beta-both LR pairs reproduced:", reproduced_pairs, "/", primary_pairs, "\n")
cat("Primary beta-supported ligands retained:", reproduced_ligands, "/", primary_ligands, "\n")
cat("Pair-level primary classification agreement:", sprintf("%.4f", pair_agreement), "\n")
cat("13-candidate beta support reproduced:", n_short_both, "/", nrow(short), "\n")
cat("13-candidate strong support reproduced:", n_short_strong_both, "/", nrow(short), "\n")

cat("\nScientific safeguards:\n")
cat("- Nearby-neuron enrichment is a spatial-region association, not neuron-intrinsic DE.\n")
cat("- Receptor compatibility is descriptive and does not demonstrate signaling.\n")
cat("- ND/T2D receptor support is not a receptor differential-expression test.\n")
cat("- The 1,600-cell dataset is a sensitivity representation, not an independent biological cohort.\n")
cat("- No ligand discovery, receptor analysis, or statistical model was rerun.\n")

cat("\nOutput figures:\n")
cat("figures/ligands/Figure_P3_DRG_Ligand_IsletCompatibility_10b.pdf\n")
cat("figures/ligands/Figure_P3_DRG_Ligand_IsletCompatibility_10b.png\n")
cat("figures/ligands/Figure_P3_DRG_Ligand_IsletCompatibility_10b.tiff\n")

cat("\n====================================================\n")
cat("SCRIPT 10b COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")
