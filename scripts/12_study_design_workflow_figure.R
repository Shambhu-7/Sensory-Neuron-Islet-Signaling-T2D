#!/usr/bin/env Rscript

# ============================================================
# Paper 3 - Figure 1
# Cross-tissue study design and analytical workflow
# Presentation-only schematic; no statistical analysis is run.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

outdir <- "figures/workflow"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_void(base_size = 11))

# Palette chosen to remain consistent with the manuscript figures.
col_drg      <- "#F8766D"
col_islet    <- "#00A9A5"
col_niche    <- "#7B2CBF"
col_trace    <- "#2C7FB8"
col_final    <- "#B3B300"
col_light    <- "#F5F5F5"
col_border   <- "#3D3D3D"
col_text     <- "#202020"
col_muted    <- "#666666"

# Canvas coordinates: x 0-16, y 0-10.
p <- ggplot() +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 10), clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(12, 18, 12, 18),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Helpers -----------------------------------------------------
box <- function(xmin, xmax, ymin, ymax, fill, border = col_border, linewidth = 0.7) {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           fill = fill, colour = border, linewidth = linewidth)
}

flow_arrow <- function(x, xend, y, yend, colour = "#777777", linewidth = 0.8) {
  annotate("segment", x = x, xend = xend, y = y, yend = yend,
           colour = colour, linewidth = linewidth,
           arrow = grid::arrow(length = unit(0.16, "inches"), type = "closed"))
}

label <- function(x, y, text, size = 3.6, fontface = "plain", colour = col_text,
                  hjust = 0.5, vjust = 0.5, lineheight = 1.08) {
  annotate("text", x = x, y = y, label = text, size = size,
           fontface = fontface, colour = colour, hjust = hjust,
           vjust = vjust, lineheight = lineheight)
}

# Title -------------------------------------------------------
p <- p +
  label(0.15, 9.70,
        "Cross-tissue strategy to prioritize candidate sensory neuron-islet signaling axes",
        size = 6.0, fontface = "bold", hjust = 0) +
  label(0.15, 9.25,
        "Human diabetic DRG spatial transcriptomics integrated with human islet single-cell data, NicheNet regulatory prediction and pancreas-traced sensory-neuron transcriptomics",
        size = 3.55, colour = col_muted, hjust = 0)

# -----------------------------------------------------------------
# TOP ROW: source datasets
# -----------------------------------------------------------------

# DRG source
p <- p +
  box(0.35, 5.05, 6.75, 8.75, "#FFF1EF", border = col_drg, linewidth = 1.0) +
  label(0.62, 8.43, "1  Human diabetic DRG", 4.25, "bold", col_drg, hjust = 0) +
  label(0.62, 8.03, "GSE295206 spatial transcriptomics", 3.55, "bold", hjust = 0) +
  label(0.62, 7.62,
        "6 donors | 7 DRGs | 16 Visium sections\n2,181 primary annotated spots",
        3.25, hjust = 0) +
  label(0.62, 7.10,
        "Nageotte-associated vs nearby-neuron regions",
        3.15, colour = col_muted, hjust = 0)

# Islet source
p <- p +
  box(5.65, 10.35, 6.75, 8.75, "#EAF9F8", border = col_islet, linewidth = 1.0) +
  label(5.92, 8.43, "2  Human pancreatic islets", 4.25, "bold", col_islet, hjust = 0) +
  label(5.92, 8.03, "GSE81608 single-cell transcriptomics", 3.55, "bold", hjust = 0) +
  label(5.92, 7.62,
        "Primary: 1,492 cells | 18 donors\nSensitivity: 1,600-cell GEO representation",
        3.25, hjust = 0) +
  label(5.92, 7.10,
        "alpha/beta receptor expression and beta-cell T2D program",
        3.15, colour = col_muted, hjust = 0)

# Trace-n-Seq source
p <- p +
  box(10.95, 15.65, 6.75, 8.75, "#EDF5FB", border = col_trace, linewidth = 1.0) +
  label(11.22, 8.43, "3  Pancreas-traced sensory neurons", 4.05, "bold", col_trace, hjust = 0) +
  label(11.22, 8.03, "E-MTAB-12940 Trace-n-Seq", 3.55, "bold", hjust = 0) +
  label(11.22, 7.62,
        "SMART-seq2 | retrograde pancreas tracing\nStrict healthy sensory cohorts",
        3.25, hjust = 0) +
  label(11.22, 7.10,
        "DRG replication plus nodose-ganglion comparison",
        3.15, colour = col_muted, hjust = 0)

# -----------------------------------------------------------------
# MIDDLE ROW: analyses
# -----------------------------------------------------------------

# DRG analysis
p <- p +
  flow_arrow(2.70, 2.70, 6.70, 6.18, colour = col_drg) +
  box(0.55, 4.85, 4.05, 6.12, "white", border = col_drg) +
  label(0.82, 5.80, "Donor-aware DRG analysis", 3.85, "bold", hjust = 0) +
  label(0.82, 5.38,
        "Donor x region pseudobulk edgeR\nCorrelation-aware redox pathway testing\nLeave-one-donor-out robustness",
        3.15, hjust = 0) +
  label(0.82, 4.52,
        "Output: degeneration-associated redox programs\nand robust curated ligand candidates",
        3.05, colour = col_muted, hjust = 0)

# Islet analysis
p <- p +
  flow_arrow(8.00, 8.00, 6.70, 6.18, colour = col_islet) +
  box(5.45, 10.55, 4.05, 6.12, "white", border = col_islet) +
  label(5.72, 5.80, "Human beta-cell receiver analysis", 3.85, "bold", hjust = 0) +
  label(5.72, 5.38,
        "High-confidence ligand-receptor compatibility\nDonor-level beta-cell T2D pseudobulk edgeR\nLODO and >=10-cell/donor sensitivity",
        3.15, hjust = 0) +
  label(5.72, 4.52,
        "Output: receptor-supported candidates\nand a 27-gene T2D receiver program",
        3.05, colour = col_muted, hjust = 0)

# Anatomical validation
p <- p +
  flow_arrow(13.30, 13.30, 6.70, 6.18, colour = col_trace) +
  box(11.15, 15.45, 4.05, 6.12, "white", border = col_trace) +
  label(11.42, 5.80, "Anatomical-expression validation", 3.75, "bold", hjust = 0) +
  label(11.42, 5.38,
        "Human-to-mouse orthologue mapping\nCandidate detection in pancreas-traced DRG\nSensory-subtype and nodose expression",
        3.15, hjust = 0) +
  label(11.42, 4.52,
        "Output: evidence that candidates are expressed\nin anatomically relevant sensory neurons",
        3.05, colour = col_muted, hjust = 0)

# -----------------------------------------------------------------
# NICHE-NET integration
# -----------------------------------------------------------------

p <- p +
  flow_arrow(4.88, 6.20, 5.08, 3.45, colour = "#888888") +
  flow_arrow(8.00, 8.00, 4.00, 3.48, colour = "#888888") +
  box(5.55, 10.45, 2.05, 3.45, "#F4ECFA", border = col_niche, linewidth = 1.0) +
  label(5.82, 3.14, "4  NicheNet regulatory compatibility", 3.95, "bold", col_niche, hjust = 0) +
  label(5.82, 2.72,
        "DRG sensory-source ligands tested against the human beta-cell T2D receiver program",
        3.05, hjust = 0) +
  label(5.82, 2.31,
        "Ligand activity = predicted regulatory compatibility, not causal signaling evidence",
        2.85, colour = col_muted, hjust = 0)

# Trace evidence feeds final integration directly
p <- p +
  flow_arrow(13.25, 11.05, 4.00, 1.58, colour = "#888888")

# NicheNet to final
p <- p +
  flow_arrow(8.00, 8.00, 2.00, 1.60, colour = col_niche)

# -----------------------------------------------------------------
# FINAL integration
# -----------------------------------------------------------------

p <- p +
  box(3.10, 12.90, 0.30, 1.58, "#FAFAE8", border = col_final, linewidth = 1.15) +
  label(3.42, 1.29, "5  Cross-tissue evidence integration", 4.15, "bold", col_final, hjust = 0) +
  label(3.42, 0.88,
        "DRG robustness + beta-cell receptor compatibility + NicheNet activity + pancreas-traced sensory-neuron expression",
        3.10, hjust = 0) +
  label(3.42, 0.49,
        "Prioritized candidates: SEMA7A   |   EDIL3   |   FLRT3   |   PGF   |   SPP1 (reference)",
        3.25, "bold", hjust = 0)

# Small framing notes
p <- p +
  label(0.35, 3.56,
        "Spatial contrast interpreted as degeneration-associated and composition-aware; not neuron-intrinsic or diabetes-induced.",
        2.72, colour = col_muted, hjust = 0) +
  label(15.65, 0.08,
        "Cross-tissue integration prioritizes testable signaling candidates; it does not demonstrate direct inter-organ communication.",
        2.55, colour = col_muted, hjust = 1)

# Save --------------------------------------------------------
png_file <- file.path(outdir, "Figure_P3_StudyDesign_Workflow_12.png")
pdf_file <- file.path(outdir, "Figure_P3_StudyDesign_Workflow_12.pdf")

# PNG
ggsave(png_file, p, width = 15.8, height = 9.4, units = "in", dpi = 300, bg = "white")

# PDF: prefer cairo for Unicode beta / en dash support when available.
if (capabilities("cairo")) {
  ggsave(pdf_file, p, width = 15.8, height = 9.4, units = "in", device = cairo_pdf, bg = "white")
} else {
  ggsave(pdf_file, p, width = 15.8, height = 9.4, units = "in", bg = "white")
}

# Source-data / figure audit record (schematic only)
audit <- data.frame(
  figure = "Figure 1",
  panel = c("Source 1", "Source 2", "Source 3", "Analysis 1", "Analysis 2", "Analysis 3", "Integration 1", "Integration 2"),
  content = c(
    "GSE295206 human diabetic DRG spatial transcriptomics: 6 donors, 7 DRGs, 16 sections, 2181 primary spots",
    "GSE81608 human islet single-cell transcriptomics: 1492-cell primary representation and 1600-cell sensitivity representation",
    "E-MTAB-12940 Trace-n-Seq pancreas-traced sensory-neuron transcriptomics",
    "Donor-aware DRG redox/pathway and ligand discovery",
    "Human beta-cell receptor compatibility and donor-level T2D receiver program",
    "Pancreas-traced sensory-neuron anatomical-expression validation",
    "NicheNet ligand activity against human beta-cell T2D receiver program",
    "Integrated prioritization: SEMA7A, EDIL3, FLRT3, PGF; SPP1 as reference signal"
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(outdir, "Figure_P3_StudyDesign_Workflow_12_Audit.csv"), row.names = FALSE)

cat("============================================================\n")
cat("SCRIPT 12 COMPLETE - FIGURE 1 STUDY DESIGN\n")
cat("============================================================\n")
cat("No statistical analysis was rerun.\n")
cat("PNG:", png_file, "\n")
cat("PDF:", pdf_file, "\n")
cat("Audit:", file.path(outdir, "Figure_P3_StudyDesign_Workflow_12_Audit.csv"), "\n")
cat("Final candidates displayed: SEMA7A, EDIL3, FLRT3, PGF; SPP1 reference.\n")
cat("============================================================\n")
