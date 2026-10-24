#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

cat("====================================================\n")
cat("SCRIPT 09b: HUMAN DIABETIC DRG REDOX PROGRAM FIGURE\n")
cat("Figure-only integration of frozen Script 02/02b outputs.\n")
cat("No edgeR, CAMERA or LODO model is rerun.\n")
cat("====================================================\n\n")

# -----------------------------
# Paths
# -----------------------------
redox_dir <- "results/redox"
fig_dir   <- "figures/redox"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

f_de <- file.path(redox_dir, "Table_P3_DE_Nageotte_vs_NearbyNeuron_edgeR.csv")
f_cam_fixed <- file.path(redox_dir, "Table_P3_CAMERA_Focused_Redox_fixedCor001.csv")
f_cam_est <- file.path(redox_dir, "Table_P3_CAMERA_Focused_Redox_estimatedCorrelation.csv")
f_lodo_path <- file.path(redox_dir, "Table_P3_LODO_CAMERA_Pathway_Robustness.csv")
f_anchor <- file.path(redox_dir, "Table_P3_Redox_Anchor_LODO_Robustness.csv")

required_files <- c(f_de, f_cam_fixed, f_cam_est, f_lodo_path, f_anchor)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required frozen input file(s):\n", paste(missing_files, collapse = "\n"))
}

# -----------------------------
# Frozen inputs
# -----------------------------
de <- read.csv(f_de, check.names = FALSE, stringsAsFactors = FALSE)
cam_fixed <- read.csv(f_cam_fixed, check.names = FALSE, stringsAsFactors = FALSE)
cam_est <- read.csv(f_cam_est, check.names = FALSE, stringsAsFactors = FALSE)
lodo_path <- read.csv(f_lodo_path, check.names = FALSE, stringsAsFactors = FALSE)
anchor <- read.csv(f_anchor, check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------
# Audit expected frozen values
# -----------------------------
stopifnot(nrow(de) == 9194)
stopifnot(sum(de$FDR < 0.05, na.rm = TRUE) == 2570)
stopifnot(sum(de$FDR < 0.05 & de$logFC > 0, na.rm = TRUE) == 1348)
stopifnot(sum(de$FDR < 0.05 & de$logFC < 0, na.rm = TRUE) == 1222)
stopifnot(sum(de$FDR < 0.05 & abs(de$logFC) >= 0.25, na.rm = TRUE) == 2553)

# -----------------------------
# Display constants
# -----------------------------
col_nageotte <- "#F8766D"
col_nearby   <- "#00BFC4"
col_ns       <- "grey82"
col_fixed    <- "#222222"
col_est      <- "grey62"

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 10, margin = margin(b = 8)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(7, 10, 7, 7)
  )

# ====================================================
# PANEL A — genome-wide differential expression
# ====================================================
de_plot <- de %>%
  mutate(
    neglog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    display_class = case_when(
      FDR < 0.05 & abs(logFC) >= 0.25 & logFC > 0 ~ "Nageotte-associated",
      FDR < 0.05 & abs(logFC) >= 0.25 & logFC < 0 ~ "Nearby-neuron-associated",
      TRUE ~ "Other"
    )
  )

# Label a compact, biologically interpretable set without changing statistics.
top_genes <- de_plot %>%
  arrange(FDR, desc(abs(logFC))) %>%
  slice_head(n = 10) %>%
  pull(gene)

redox_label_genes <- c("NFE2L2", "TREM2", "TYROBP", "TXN", "PRDX2", "SRXN1", "ACSL4", "SPP1")
label_genes <- unique(c(top_genes, redox_label_genes))
label_df <- de_plot %>% filter(gene %in% label_genes)

pA <- ggplot(de_plot, aes(x = logFC, y = neglog10FDR)) +
  geom_point(
    aes(color = display_class),
    size = 1.15, alpha = 0.62, stroke = 0
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", linewidth = 0.35, color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35, color = "grey55") +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.0,
    min.segment.length = 0,
    box.padding = 0.28,
    point.padding = 0.18,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 295206,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Nageotte-associated" = col_nageotte,
      "Nearby-neuron-associated" = col_nearby,
      "Other" = col_ns
    ),
    breaks = c("Nageotte-associated", "Nearby-neuron-associated", "Other")
  ) +
  labs(
    title = "A  Degeneration-associated transcriptional landscape",
    subtitle = "Donor-adjusted pseudobulk edgeR; 6 donors, 2,181 primary Visium spots",
    x = expression(log[2]~fold~change~"(Nageotte - nearby neuron)"),
    y = expression(-log[10]~FDR),
    color = NULL
  ) +
  base_theme +
  theme(legend.position = "bottom")

# ====================================================
# PANEL B — CAMERA correlation sensitivity
# ====================================================
pathways <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "REACTOME_MITOPHAGY",
  "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2",
  "GOBP_NEUROINFLAMMATORY_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)

path_labels <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNF-alpha / NF-kB signaling",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammatory response",
  "REACTOME_MITOPHAGY" = "Mitophagy",
  "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2" = "NFE2L2 nuclear events",
  "GOBP_NEUROINFLAMMATORY_RESPONSE" = "Neuroinflammatory response",
  "HALLMARK_APOPTOSIS" = "Apoptosis",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" = "ROS pathway"
)

cam_b <- bind_rows(
  cam_fixed %>% mutate(Model = "Fixed correlation (0.01)"),
  cam_est %>% mutate(Model = "Estimated correlation")
) %>%
  filter(gene_set %in% pathways) %>%
  mutate(
    pathway = unname(path_labels[gene_set]),
    pathway = factor(pathway, levels = rev(unname(path_labels[pathways]))),
    signed_score = ifelse(Direction == "Up", 1, -1) * -log10(pmax(FDR, .Machine$double.xmin)),
    Model = factor(Model, levels = c("Fixed correlation (0.01)", "Estimated correlation"))
  )

# Connect the two sensitivity estimates for each pathway.
pB <- ggplot(cam_b, aes(x = signed_score, y = pathway)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey35") +
  geom_vline(xintercept = c(-(-log10(0.05)), -log10(0.05)),
             linetype = "dashed", linewidth = 0.35, color = "grey60") +
  geom_line(aes(group = pathway), color = "grey70", linewidth = 0.6) +
  geom_point(aes(shape = Model, color = Direction), size = 3.3, stroke = 0.9) +
  scale_shape_manual(
    values = c("Fixed correlation (0.01)" = 16, "Estimated correlation" = 17),
    guide = guide_legend(override.aes = list(color = "black"))
  ) +
  scale_color_manual(
    values = c("Up" = col_nageotte, "Down" = col_nearby),
    breaks = c("Down", "Up"),
    labels = c("Down" = "Nearby-neuron-associated", "Up" = "Nageotte-associated"),
    guide = guide_legend(override.aes = list(shape = 16, size = 3.3))
  ) +
  labs(
    title = "B  Redox and inflammatory pathway programs",
    subtitle = "CAMERA sensitivity to inter-gene correlation assumptions",
    x = expression("Signed " * -log[10] * " FDR"),
    y = NULL,
    shape = "CAMERA model",
    color = "Direction"
  ) +
  base_theme +
  theme(
    axis.text.y = element_text(size = 9.5),
    legend.position = "bottom",
    legend.box = "vertical"
  )

# ====================================================
# PANEL C — leave-one-donor-out pathway robustness
# ====================================================
lodo_c <- lodo_path %>%
  filter(gene_set %in% pathways) %>%
  transmute(
    gene_set,
    pathway = unname(path_labels[gene_set]),
    `Fixed: direction` = fixed_same_direction,
    `Fixed: FDR < 0.05` = fixed_FDR05,
    `Estimated: direction` = estimated_same_direction,
    `Estimated: FDR < 0.05` = estimated_FDR05
  ) %>%
  pivot_longer(
    cols = c(`Fixed: direction`, `Fixed: FDR < 0.05`,
             `Estimated: direction`, `Estimated: FDR < 0.05`),
    names_to = "metric", values_to = "n_of_6"
  ) %>%
  mutate(
    pathway = factor(pathway, levels = rev(unname(path_labels[pathways]))),
    metric = factor(
      metric,
      levels = c("Fixed: direction", "Fixed: FDR < 0.05",
                 "Estimated: direction", "Estimated: FDR < 0.05")
    ),
    text_col = ifelse(n_of_6 <= 2, "white", "black")
  )

pC <- ggplot(lodo_c, aes(x = metric, y = pathway, fill = n_of_6)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = n_of_6, color = text_col), size = 3.4, fontface = "bold") +
  scale_color_identity() +
  scale_fill_viridis_c(limits = c(0, 6), breaks = 0:6) +
  labs(
    title = "C  Donor robustness of pathway direction",
    subtitle = "Number of six leave-one-donor-out fits meeting each criterion",
    x = NULL,
    y = NULL,
    fill = "LODO fits"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1, size = 8.8),
    axis.text.y = element_text(size = 9.5),
    legend.position = "bottom"
  )

# ====================================================
# PANEL D — selected redox / inflammatory anchor genes
# ====================================================
anchor_genes <- c(
  "NFE2L2", "FTH1", "FTL", "TYROBP", "TREM2",
  "GSS", "TXN", "PRDX2", "SRXN1", "SOD1", "ACSL4", "CLU", "SPP1"
)

anchor_d <- anchor %>%
  filter(gene %in% anchor_genes) %>%
  mutate(
    gene = factor(gene, levels = rev(anchor_genes)),
    region = ifelse(logFC > 0, "Nageotte-associated", "Nearby-neuron-associated"),
    neglog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    donor_concordance = ifelse(all_6_donors_same_direction, "6/6 donor direction", "<6/6 donor direction"),
    lodo_label = paste0(n_LODO_FDR05, "/6")
  )

if (nrow(anchor_d) != length(anchor_genes)) {
  stop("Not all prespecified anchor genes were found in the frozen anchor table.")
}

pD <- ggplot(anchor_d, aes(x = logFC, y = gene)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey35") +
  geom_segment(aes(x = 0, xend = logFC, yend = gene, color = region), linewidth = 0.75) +
  geom_point(
    aes(color = region, size = neglog10FDR, shape = donor_concordance),
    stroke = 0.65
  ) +
  geom_text(
    aes(label = lodo_label),
    hjust = ifelse(anchor_d$logFC >= 0, -0.35, 1.35),
    size = 2.8,
    color = "grey25"
  ) +
  scale_color_manual(values = c(
    "Nageotte-associated" = col_nageotte,
    "Nearby-neuron-associated" = col_nearby
  )) +
  scale_shape_manual(values = c("6/6 donor direction" = 16, "<6/6 donor direction" = 1)) +
  scale_size_continuous(range = c(2.5, 6.0)) +
  scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
  labs(
    title = "D  Selected redox and degeneration-associated genes",
    subtitle = "Labels show LODO fits with FDR < 0.05 (of 6); shape shows donor-direction concordance",
    x = expression(log[2]~fold~change~"(Nageotte - nearby neuron)"),
    y = NULL,
    color = "Region",
    size = expression(-log[10]~FDR),
    shape = NULL
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.y = element_text(face = "italic")
  )

# -----------------------------
# Source-data freeze for Figure 2
# -----------------------------
write.csv(
  de_plot %>% select(gene, logFC, logCPM, F, PValue, FDR, direction, display_class),
  file.path(redox_dir, "Figure09_SourceData_PanelA_Volcano.csv"),
  row.names = FALSE
)

write.csv(
  cam_b %>% select(gene_set, pathway, NGenes, Direction, PValue, FDR, Model, signed_score),
  file.path(redox_dir, "Figure09_SourceData_PanelB_CAMERA.csv"),
  row.names = FALSE
)

write.csv(
  lodo_c %>% select(gene_set, pathway, metric, n_of_6),
  file.path(redox_dir, "Figure09_SourceData_PanelC_LODO.csv"),
  row.names = FALSE
)

write.csv(
  anchor_d %>% select(gene, logFC, FDR, direction, n_LODO_tested,
                      n_LODO_same_direction, n_LODO_FDR05,
                      all_LODO_direction_concordant,
                      n_donors_same_direction, all_6_donors_same_direction),
  file.path(redox_dir, "Figure09_SourceData_PanelD_Anchors.csv"),
  row.names = FALSE
)

# -----------------------------
# Export helpers
# -----------------------------
save_panel <- function(plot_obj, stem, width = 7.0, height = 5.2) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot_obj,
         device = cairo_pdf, width = width, height = height, units = "in")
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot_obj,
         width = width, height = height, units = "in", dpi = 300, bg = "white")
}

save_panel(pA, "Figure09_A_Volcano", 7.0, 5.4)
save_panel(pB, "Figure09_B_CAMERA", 7.0, 5.4)
save_panel(pC, "Figure09_C_LODO_Pathways", 7.0, 5.4)
save_panel(pD, "Figure09_D_RedoxAnchors", 7.0, 5.6)

# -----------------------------
# Combined Figure 2
# -----------------------------
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Degeneration-associated redox and inflammatory programs in human diabetic DRG",
    subtitle = "Donor-aware spatial pseudobulk analysis with correlation-aware pathway testing and leave-one-donor-out robustness",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      plot.subtitle = element_text(size = 11.5, hjust = 0, margin = margin(b = 8)),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

ggsave(
  file.path(fig_dir, "Figure_P3_DRG_RedoxPrograms_09b.pdf"),
  combined,
  device = cairo_pdf,
  width = 14.0, height = 10.6, units = "in"
)

ggsave(
  file.path(fig_dir, "Figure_P3_DRG_RedoxPrograms_09b.png"),
  combined,
  width = 14.0, height = 10.6, units = "in", dpi = 300, bg = "white"
)

ggsave(
  file.path(fig_dir, "Figure_P3_DRG_RedoxPrograms_09b.tiff"),
  combined,
  width = 14.0, height = 10.6, units = "in",
  dpi = 600, compression = "lzw", bg = "white"
)

# -----------------------------
# Console audit
# -----------------------------
cat("\n====================================================\n")
cat("FIGURE 09 AUDIT\n")
cat("====================================================\n")
cat("Genes tested:", nrow(de), "\n")
cat("DE FDR < 0.05:", sum(de$FDR < 0.05, na.rm = TRUE), "\n")
cat("Nageotte-associated DE:", sum(de$FDR < 0.05 & de$logFC > 0, na.rm = TRUE), "\n")
cat("Nearby-neuron-associated DE:", sum(de$FDR < 0.05 & de$logFC < 0, na.rm = TRUE), "\n")
cat("FDR < 0.05 and |logFC| >= 0.25:",
    sum(de$FDR < 0.05 & abs(de$logFC) >= 0.25, na.rm = TRUE), "\n")
cat("Panel B pathways:", length(unique(cam_b$gene_set)), "\n")
cat("Panel C pathways:", length(unique(lodo_c$gene_set)), "\n")
cat("Panel D anchors:", nrow(anchor_d), "\n")

cat("\nFixed-correlation FDR < 0.05 among displayed pathways:\n")
print(
  cam_b %>%
    filter(Model == "Fixed correlation (0.01)", FDR < 0.05) %>%
    select(pathway, Direction, FDR),
  row.names = FALSE
)

cat("\nEstimated-correlation directions/FDR among displayed pathways:\n")
print(
  cam_b %>%
    filter(Model == "Estimated correlation") %>%
    select(pathway, Direction, FDR),
  row.names = FALSE
)

cat("\nPanel D anchors:\n")
print(
  anchor_d %>%
    arrange(gene) %>%
    select(gene, logFC, FDR, n_LODO_FDR05,
           n_donors_same_direction, all_6_donors_same_direction),
  row.names = FALSE
)

cat("\nScientific safeguards:\n")
cat("- Contrast is Nageotte-associated versus nearby-neuron-associated spatial regions.\n")
cat("- Positive logFC denotes Nageotte-associated enrichment.\n")
cat("- Negative logFC denotes nearby-neuron-associated enrichment.\n")
cat("- Spatial-region differences are composition-confounded and are not called neuron-intrinsic.\n")
cat("- No normal human DRG control exists in GSE295206; results are not called diabetes-induced.\n")
cat("- Fixed-correlation CAMERA is shown together with estimated-correlation sensitivity.\n")
cat("- LODO values are robustness summaries from frozen Script 02b outputs.\n")
cat("- No new statistical analysis was performed by Script 09.\n")

cat("\nOutput figures:\n")
cat(file.path(fig_dir, "Figure_P3_DRG_RedoxPrograms_09b.pdf"), "\n")
cat(file.path(fig_dir, "Figure_P3_DRG_RedoxPrograms_09b.png"), "\n")
cat(file.path(fig_dir, "Figure_P3_DRG_RedoxPrograms_09b.tiff"), "\n")

cat("\n====================================================\n")
cat("SCRIPT 09b COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")
