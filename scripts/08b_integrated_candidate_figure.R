suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

cat("====================================================\n")
cat("SCRIPT 08b: PUBLICATION-POLISHED INTEGRATED FIGURE\n")
cat("Presentation-only revision; no analytical values are changed.\n")
cat("====================================================\n\n")

dir.create("figures/integration", recursive = TRUE, showWarnings = FALSE)
dir.create("results/integration/figure_source_data", recursive = TRUE, showWarnings = FALSE)
dir.create("reproducibility", recursive = TRUE, showWarnings = FALSE)

f_main <- "results/integration/Table_P3_Final_MainCandidate_Hierarchy.csv"
f_sub  <- "results/integration/Table_P3_Final_SensorySubtypeEvidence.csv"

if (!file.exists(f_main)) stop("Missing: ", f_main)
if (!file.exists(f_sub)) stop("Missing: ", f_sub)

x <- read.csv(f_main, check.names = FALSE, stringsAsFactors = FALSE)
sub <- read.csv(f_sub, check.names = FALSE, stringsAsFactors = FALSE)

genes <- c("SEMA7A", "EDIL3", "FLRT3", "PGF", "SPP1")
if (!all(genes %in% x$HumanGene)) stop("Frozen five-candidate hierarchy is incomplete.")

x <- x[match(genes, x$HumanGene), , drop = FALSE]
lev <- rev(genes)
x$HumanGene <- factor(x$HumanGene, levels = lev)

theme_paper <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9.2),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 10, 8, 8)
  )

# ------------------------------------------------------------
# PANEL A — human diabetic DRG spatial evidence
# Negative original logFC = nearby-neuron enrichment.
# Display -logFC so positive values indicate nearby-neuron enrichment.
# ------------------------------------------------------------

A <- data.frame(
  HumanGene = x$HumanGene,
  NearbyNeuron_Enrichment = -x$LR_logFC,
  FDR = x$LR_FDR,
  Tier = x$LR_candidate_tier,
  stringsAsFactors = FALSE
)

A$MinusLog10FDR <- -log10(pmax(A$FDR, .Machine$double.xmin))
A$Tier <- ifelse(
  A$Tier == "Tier_A_strict", "Tier A (strict)",
  ifelse(A$Tier == "Tier_B_LODO_robust", "Tier B (LODO-robust)", A$Tier)
)
A$Tier <- factor(A$Tier, levels = c("Tier A (strict)", "Tier B (LODO-robust)"))

write.csv(
  A,
  "results/integration/figure_source_data/Figure_Integrated_08b_PanelA_DRG.csv",
  row.names = FALSE
)

pA <- ggplot(A, aes(x = NearbyNeuron_Enrichment, y = HumanGene)) +
  geom_segment(
    aes(x = 0, xend = NearbyNeuron_Enrichment, yend = HumanGene),
    linewidth = 0.65
  ) +
  geom_point(aes(size = MinusLog10FDR, shape = Tier)) +
  scale_size_continuous(
    name = expression(-log[10]~FDR),
    range = c(3.0, 6.6)
  ) +
  labs(
    title = "A  Human diabetic DRG",
    subtitle = "Nearby-neuron-associated ligand enrichment",
    x = expression("Nearby-neuron enrichment (" * -log[2] * "FC)"),
    y = NULL,
    shape = "DRG tier"
  ) +
  coord_cartesian(
    xlim = c(0, max(A$NearbyNeuron_Enrichment, na.rm = TRUE) * 1.12),
    clip = "off"
  ) +
  theme_paper +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 1),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    shape = guide_legend(order = 1, nrow = 1),
    size = guide_legend(order = 2, nrow = 1)
  )

# ------------------------------------------------------------
# PANEL B — human beta-cell receptor compatibility
# ------------------------------------------------------------

B <- rbind(
  data.frame(
    HumanGene = as.character(x$HumanGene),
    Support = "Supported in ND + T2D",
    N_receptors = x$LR_n_HC_beta_both_receptors,
    stringsAsFactors = FALSE
  ),
  data.frame(
    HumanGene = as.character(x$HumanGene),
    Support = "Strong in ND + T2D",
    N_receptors = x$LR_n_HC_beta_strong_both_receptors,
    stringsAsFactors = FALSE
  )
)

B$HumanGene <- factor(B$HumanGene, levels = lev)
B$Support <- factor(
  B$Support,
  levels = c("Supported in ND + T2D", "Strong in ND + T2D")
)

write.csv(
  B,
  "results/integration/figure_source_data/Figure_Integrated_08b_PanelB_Receptors.csv",
  row.names = FALSE
)

pB <- ggplot(B, aes(x = N_receptors, y = HumanGene, fill = Support)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_x_continuous(
    breaks = 0:5,
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "B  Human \u03b2-cell receptors",
    subtitle = "High-confidence receptor compatibility",
    x = "Number of \u03b2-cell receptors",
    y = NULL,
    fill = NULL
  ) +
  theme_paper +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 1)
  ) +
  guides(fill = guide_legend(nrow = 1))

# ------------------------------------------------------------
# PANEL C — NicheNet receiver-program prediction
# ------------------------------------------------------------

C <- rbind(
  data.frame(
    HumanGene = as.character(x$HumanGene),
    Receiver = "Primary (22/27)",
    Corrected_AUPR = x$NicheNet_Primary_AUPR_corrected,
    Rank = x$NicheNet_Primary_Rank,
    stringsAsFactors = FALSE
  ),
  data.frame(
    HumanGene = as.character(x$HumanGene),
    Receiver = "Core (8/10)",
    Corrected_AUPR = x$NicheNet_Core_AUPR_corrected,
    Rank = x$NicheNet_Core_Rank,
    stringsAsFactors = FALSE
  )
)

C$HumanGene <- factor(C$HumanGene, levels = lev)
C$Receiver <- factor(C$Receiver, levels = c("Primary (22/27)", "Core (8/10)"))

write.csv(
  C,
  "results/integration/figure_source_data/Figure_Integrated_08b_PanelC_NicheNet.csv",
  row.names = FALSE
)

rank_nudge <- max(C$Corrected_AUPR, na.rm = TRUE) * 0.035

pC <- ggplot(C, aes(x = Corrected_AUPR, y = HumanGene)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(aes(group = HumanGene), linewidth = 0.6) +
  geom_point(aes(shape = Receiver), size = 3.1) +
  geom_text(
    aes(label = paste0("#", Rank)),
    nudge_x = rank_nudge,
    hjust = 0,
    size = 3.0,
    show.legend = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.20))) +
  labs(
    title = "C  \u03b2-cell receiver prediction",
    subtitle = "NicheNet ligand activity",
    x = "Corrected AUPR",
    y = NULL,
    shape = NULL
  ) +
  theme_paper +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 1)
  ) +
  guides(shape = guide_legend(nrow = 1))

# ------------------------------------------------------------
# PANEL D — Trace-n-Seq anatomical validation
# ------------------------------------------------------------

D <- rbind(
  data.frame(
    HumanGene = as.character(x$HumanGene),
    Cohort = "DRG batch 1",
    PctDetected = x$DRG_Batch1_PctDetected,
    stringsAsFactors = FALSE
  ),
  data.frame(
    HumanGene = as.character(x$HumanGene),
    Cohort = "DRG batch 2",
    PctDetected = x$DRG_Batch2_PctDetected,
    stringsAsFactors = FALSE
  ),
  data.frame(
    HumanGene = as.character(x$HumanGene),
    Cohort = "NG",
    PctDetected = x$NG_PctDetected,
    stringsAsFactors = FALSE
  )
)

D$HumanGene <- factor(D$HumanGene, levels = lev)
D$Cohort <- factor(D$Cohort, levels = c("DRG batch 1", "DRG batch 2", "NG"))
D$TextClass <- ifelse(D$PctDetected < 35, "light", "dark")

write.csv(
  D,
  "results/integration/figure_source_data/Figure_Integrated_08b_PanelD_TraceNSeq.csv",
  row.names = FALSE
)

pD <- ggplot(D, aes(x = Cohort, y = HumanGene, fill = PctDetected)) +
  geom_tile(linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.0f", PctDetected), color = TextClass),
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("light" = "white", "dark" = "black")) +
  scale_fill_viridis_c(limits = c(0, 100), name = "% detected") +
  labs(
    title = "D  Pancreas-traced neurons",
    subtitle = "Trace-n-Seq anatomical validation",
    x = NULL,
    y = NULL
  ) +
  theme_paper +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(
    title.position = "left",
    title.hjust = 1,
    barwidth = unit(1.45, "in"),
    barheight = unit(0.13, "in")
  ))

# ------------------------------------------------------------
# PANEL E — pancreas-traced DRG sensory subtypes
# Descriptive detection percentages; not formal enrichment testing.
# ------------------------------------------------------------

need <- c("HumanGene", "CellType", "N_cells", "Pct_detected", "Mean_log1p_CPM")
if (!all(need %in% names(sub))) {
  stop(
    "Subtype table lacks expected columns: ",
    paste(setdiff(need, names(sub)), collapse = ", ")
  )
}

sub <- sub[sub$HumanGene %in% genes, , drop = FALSE]

short <- c(
  "Peripheral sensory neurofilament neurons" = "Neurofilament",
  "Peripheral sensory non-peptidergic neurons" = "Non-peptidergic",
  "Peripheral sensory peptidergic neurons" = "Peptidergic"
)

sub$SensorySubtype <- unname(short[as.character(sub$CellType)])
sub <- sub[!is.na(sub$SensorySubtype), , drop = FALSE]

sub$HumanGene <- factor(sub$HumanGene, levels = lev)
sub$SensorySubtype <- factor(
  sub$SensorySubtype,
  levels = c("Neurofilament", "Non-peptidergic", "Peptidergic")
)
sub$TextClass <- ifelse(sub$Pct_detected < 35, "light", "dark")

write.csv(
  sub,
  "results/integration/figure_source_data/Figure_Integrated_08b_PanelE_Subtypes.csv",
  row.names = FALSE
)

pE <- ggplot(sub, aes(x = SensorySubtype, y = HumanGene, fill = Pct_detected)) +
  geom_tile(linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.0f", Pct_detected), color = TextClass),
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("light" = "white", "dark" = "black")) +
  scale_fill_viridis_c(limits = c(0, 100), name = "% detected") +
  labs(
    title = "E  Sensory-neuron subtype",
    subtitle = "Combined pancreas-traced DRG cohorts",
    x = NULL,
    y = NULL
  ) +
  theme_paper +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(
    title.position = "left",
    title.hjust = 1,
    barwidth = unit(1.45, "in"),
    barheight = unit(0.13, "in")
  ))

# ------------------------------------------------------------
# COMBINE
# ------------------------------------------------------------

top <- pA + pB + pC +
  plot_layout(widths = c(1.10, 0.95, 1.10))

bottom <- pD + pE +
  plot_layout(widths = c(1, 1))

combined <- top / bottom +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(
    title = "Integrated prioritization of candidate sensory neuron\u2013islet signaling axes",
    subtitle = paste0(
      "Cross-tissue evidence from human diabetic DRG, human \u03b2 cells, ",
      "NicheNet and pancreas-traced sensory neurons"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9.5),
      plot.margin = margin(8, 10, 6, 10)
    )
  )

# ------------------------------------------------------------
# EXPORT INDIVIDUAL POLISHED PANELS
# ------------------------------------------------------------

panels <- list(
  A_DRG = pA,
  B_Receptors = pB,
  C_NicheNet = pC,
  D_TraceNSeq = pD,
  E_Subtypes = pE
)

for (nm in names(panels)) {
  ggsave(
    paste0("figures/integration/Figure_Integrated_08b_", nm, ".pdf"),
    panels[[nm]],
    device = cairo_pdf,
    width = 6.4, height = 4.4, units = "in"
  )
  ggsave(
    paste0("figures/integration/Figure_Integrated_08b_", nm, ".png"),
    panels[[nm]],
    width = 6.4, height = 4.4, units = "in", dpi = 400
  )
}

# ------------------------------------------------------------
# EXPORT COMBINED POLISHED FIGURE
# ------------------------------------------------------------

ggsave(
  "figures/integration/Figure_P3_IntegratedCandidateEvidence_08b.pdf",
  combined,
  device = cairo_pdf,
  width = 13.5, height = 9.2, units = "in"
)

ggsave(
  "figures/integration/Figure_P3_IntegratedCandidateEvidence_08b.png",
  combined,
  width = 13.5, height = 9.2, units = "in", dpi = 400
)

ggsave(
  "figures/integration/Figure_P3_IntegratedCandidateEvidence_08b.tiff",
  combined,
  width = 13.5, height = 9.2, units = "in", dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# AUDIT: verify no analytical values changed
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("08b FIGURE AUDIT\n")
cat("====================================================\n")
cat("Panel A candidates:", nrow(A), "\n")
cat("Panel B rows:", nrow(B), "\n")
cat("Panel C rows:", nrow(C), "\n")
cat("Panel D rows:", nrow(D), "\n")
cat("Panel E rows:", nrow(sub), "\n")
cat("Candidate order:", paste(genes, collapse = " -> "), "\n")
cat("Main-table rows used:", nrow(x), "\n")

cat("\nPresentation changes only:\n")
cat("- publication-friendly Tier labels\n")
cat("- beta replaced by \u03b2 in display labels\n")
cat("- shortened NicheNet legend labels\n")
cat("- compact legends and adjusted panel spacing\n")
cat("- contrast-aware heatmap text\n")
cat("- shortened overall subtitle\n")

cat("\nScientific safeguards retained:\n")
cat("- DRG contrast remains nearby-neuron-associated, not diabetes-induced\n")
cat("- receptor compatibility remains descriptive\n")
cat("- Trace-n-Seq remains anatomical-expression validation\n")
cat("- NG remains labelled NG\n")
cat("- no composite score or new candidate ranking\n")
cat("- no statistical analysis rerun\n")

writeLines(
  capture.output(sessionInfo()),
  "reproducibility/08b_integrated_candidate_figure_sessionInfo.txt"
)

cat("\nOutput figures:\n")
cat("figures/integration/Figure_P3_IntegratedCandidateEvidence_08b.pdf\n")
cat("figures/integration/Figure_P3_IntegratedCandidateEvidence_08b.tiff\n")
cat("figures/integration/Figure_P3_IntegratedCandidateEvidence_08b.png\n")

cat("\n====================================================\n")
cat("SCRIPT 08b COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")
