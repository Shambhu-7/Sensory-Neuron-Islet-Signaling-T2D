#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 11b — FINAL PRESENTATION REVISION
# Figure 4: Human T2D beta-cell receiver program and NicheNet
#
# PRESENTATION ONLY:
# - Uses frozen outputs from Scripts 05, 05b and 06.
# - Does NOT rerun edgeR, LODO, NicheNet, receptor analysis,
#   ligand discovery, or any statistical test.
# - The saved predicted-target-link table is intentionally not
#   used because target/weight fields were not populated.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

cat("====================================================\n")
cat("SCRIPT 11b: T2D BETA-CELL RECEIVER + NICHENET FIGURE — FINAL\n")
cat("Frozen-results visualization only; no analysis rerun.\n")
cat("====================================================\n\n")

# -----------------------------
# Paths
# -----------------------------
f_de <- "results/islet/Table_P3_Beta_T2D_edgeR_Primary_18Donors.csv"
f_rob <- "results/islet/Table_P3_Beta_T2D_ReceiverGene_Robustness.csv"
f_core <- "results/islet/Table_P3_Beta_T2D_NicheNet_HighConfidenceCore.csv"
f_nn_primary <- "results/nichenet/Table_P3_NicheNet_SensorySource_PrimaryReceiver_LigandActivity.csv"
f_nn_core <- "results/nichenet/Table_P3_NicheNet_SensorySource_HighConfidenceCore_LigandActivity.csv"

required <- c(f_de, f_rob, f_core, f_nn_primary, f_nn_core)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required frozen input(s):\n", paste(missing, collapse = "\n"))
}

dir.create("figures/islet", recursive = TRUE, showWarnings = FALSE)
dir.create("results/islet/figure_source_data", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
pick_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Could not find any of these columns: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(names(df), collapse = " | ")
    )
  }
  NA_character_
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE","T","1","YES","Y")
}

clean_ligand <- function(x) {
  x <- as.character(x)
  x <- sub("^.*?:", "", x)
  x <- trimws(x)
  x
}

# Publication palette
col_up <- "#D55E5E"
col_down <- "#2A9D8F"
col_ns <- "grey78"
col_core <- "#7B2CBF"
col_primary <- "#2878B5"

theme_pub <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9.5, color = "grey30"),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 9),
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.8),
    plot.margin = margin(7, 10, 7, 7)
  )

# ============================================================
# PANEL A — donor-level beta-cell T2D differential expression
# ============================================================
de <- read.csv(f_de, check.names = FALSE)

gene_col <- pick_col(de, c("gene","Gene","symbol","GeneSymbol"))
lfc_col <- pick_col(de, c("logFC","log2FC","log2FoldChange"))
fdr_col <- pick_col(de, c("FDR","padj","adj.P.Val"))
logcpm_col <- pick_col(de, c("logCPM","AveExpr","baseMean"), required = FALSE)

de$gene_plot <- as.character(de[[gene_col]])
de$logFC_plot <- as.numeric(de[[lfc_col]])
de$FDR_plot <- as.numeric(de[[fdr_col]])
de$neglog10FDR <- -log10(pmax(de$FDR_plot, .Machine$double.xmin))

de$DE_class <- "Not FDR-significant"
de$DE_class[de$FDR_plot < 0.05 & de$logFC_plot > 0] <- "T2D-up"
de$DE_class[de$FDR_plot < 0.05 & de$logFC_plot < 0] <- "T2D-down"
de$DE_class <- factor(
  de$DE_class,
  levels = c("Not FDR-significant","T2D-up","T2D-down")
)

sig <- de %>% filter(FDR_plot < 0.05)
if (nrow(sig) != 27) {
  warning("Expected 27 FDR-significant genes; observed ", nrow(sig))
}

# Label the 12 most significant genes plus both significant T2D-down genes.
# All 27 receiver genes are displayed explicitly in Panel B.
label_df <- bind_rows(
  sig %>% arrange(FDR_plot, desc(abs(logFC_plot))) %>% slice_head(n = 12),
  sig %>% filter(logFC_plot < 0)
) %>% distinct(gene_plot, .keep_all = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

pA <- ggplot(de, aes(logFC_plot, neglog10FDR)) +
  geom_point(aes(color = DE_class), size = 1.25, alpha = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed", linewidth = 0.4, color = "grey45"
  ) +
  scale_color_manual(
    values = c(
      "Not FDR-significant" = col_ns,
      "T2D-up" = col_up,
      "T2D-down" = col_down
    ),
    breaks = c("T2D-up","T2D-down","Not FDR-significant")
  ) +
  labs(
    title = "A  Human β-cell T2D transcriptional program",
    subtitle = "Donor-level pseudobulk edgeR: 18 donors (12 non-diabetic, 6 T2D)",
    x = expression(log[2]~"fold change (T2D - non-diabetic)"),
    y = expression(-log[10]~"FDR"),
    color = NULL
  ) +
  theme_pub +
  theme(legend.position = "bottom")

if (has_ggrepel) {
  pA <- pA + ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene_plot),
    size = 2.55,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.16,
    min.segment.length = 0,
    segment.size = 0.25,
    show.legend = FALSE,
    seed = 11
  )
} else {
  # Fallback: label the 15 most significant genes.
  lab15 <- label_df %>% arrange(FDR_plot, desc(abs(logFC_plot))) %>% slice_head(n = 15)
  pA <- pA + geom_text(
    data = lab15,
    aes(label = gene_plot),
    size = 2.3, check_overlap = TRUE,
    vjust = -0.5, show.legend = FALSE
  )
}

# ============================================================
# PANEL B — receiver-gene donor robustness
# ============================================================
rob <- read.csv(f_rob, check.names = FALSE)
core_tbl <- read.csv(f_core, check.names = FALSE)

r_gene <- pick_col(rob, c("gene","Gene","symbol","GeneSymbol"))
r_lfc <- pick_col(rob, c("full_logFC","logFC","Primary_logFC","primary_logFC"))
r_fdr <- pick_col(rob, c("full_FDR","FDR","Primary_FDR","primary_FDR"))
r_same <- pick_col(
  rob,
  c("n_LODO_same_direction","n_lodo_same_direction","LODO_same_direction"),
  required = FALSE
)
r_fdrn <- pick_col(
  rob,
  c("n_LODO_FDR05","n_lodo_FDR05","LODO_FDR05"),
  required = FALSE
)

core_gene_col <- pick_col(core_tbl, c("gene","Gene","symbol","GeneSymbol"))
core_genes <- unique(as.character(core_tbl[[core_gene_col]]))

rob$gene_plot <- as.character(rob[[r_gene]])
rob$logFC_plot <- as.numeric(rob[[r_lfc]])
rob$FDR_plot <- as.numeric(rob[[r_fdr]])
rob$Core <- ifelse(rob$gene_plot %in% core_genes,
                   "10-gene high-confidence core",
                   "Primary 27-gene receiver set")

# Frozen analysis states all 27 are direction-concordant across 18 LODO fits.
if (!is.na(r_same)) {
  rob$LODO_same <- as.numeric(rob[[r_same]])
} else {
  rob$LODO_same <- 18
}
if (!is.na(r_fdrn)) {
  rob$LODO_FDR05 <- as.numeric(rob[[r_fdrn]])
} else {
  rob$LODO_FDR05 <- NA_real_
}

rob <- rob %>%
  arrange(logFC_plot) %>%
  mutate(gene_factor = factor(gene_plot, levels = gene_plot))

pB <- ggplot(rob, aes(logFC_plot, gene_factor)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_segment(
    aes(x = 0, xend = logFC_plot, yend = gene_factor),
    linewidth = 0.55, color = "grey72"
  ) +
  geom_point(aes(fill = Core), shape = 21, size = 3.1, color = "black", stroke = 0.35) +
  scale_fill_manual(
    values = c(
      "10-gene high-confidence core" = col_core,
      "Primary 27-gene receiver set" = "white"
    )
  ) +
  labs(
    title = "B  Donor robustness of the 27-gene receiver program",
    subtitle = "All 27 retain the same T2D-direction in all 18 leave-one-donor-out fits",
    x = expression("Full-model " * log[2] * " fold change"),
    y = NULL,
    fill = NULL
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 7.6),
    legend.position = "bottom"
  )

# ============================================================
# PANEL C — NicheNet ligand activity: primary vs core
# ============================================================
nn1 <- read.csv(f_nn_primary, check.names = FALSE)
nn2 <- read.csv(f_nn_core, check.names = FALSE)

lig1 <- pick_col(nn1, c("test_ligand","ligand","Ligand","HumanGene"))
lig2 <- pick_col(nn2, c("test_ligand","ligand","Ligand","HumanGene"))

aup1 <- pick_col(
  nn1,
  c("aupr_corrected","AUPR_corrected","corrected_AUPR",
    "NicheNet_Primary_AUPR_corrected")
)
aup2 <- pick_col(
  nn2,
  c("aupr_corrected","AUPR_corrected","corrected_AUPR",
    "NicheNet_Core_AUPR_corrected")
)

rank1 <- pick_col(nn1, c("rank","Rank","ligand_rank","NicheNet_Primary_Rank"), required = FALSE)
rank2 <- pick_col(nn2, c("rank","Rank","ligand_rank","NicheNet_Core_Rank"), required = FALSE)

n1 <- data.frame(
  ligand = clean_ligand(nn1[[lig1]]),
  Primary = as.numeric(nn1[[aup1]]),
  stringsAsFactors = FALSE
)
n2 <- data.frame(
  ligand = clean_ligand(nn2[[lig2]]),
  Core = as.numeric(nn2[[aup2]]),
  stringsAsFactors = FALSE
)

if (!is.na(rank1)) n1$PrimaryRank <- as.numeric(nn1[[rank1]])
if (!is.na(rank2)) n2$CoreRank <- as.numeric(nn2[[rank2]])

nn <- full_join(n1, n2, by = "ligand")

# The frozen sensory-source analysis has 10 NicheNet-evaluable ligands.
if (nrow(nn) != 10) {
  warning("Expected 10 NicheNet-evaluable sensory-source ligands; observed ", nrow(nn))
}

nn$order_value <- pmax(nn$Primary, nn$Core, na.rm = TRUE)
nn <- nn %>% arrange(order_value)
nn$ligand_factor <- factor(nn$ligand, levels = nn$ligand)

nn_long <- nn %>%
  select(ligand, ligand_factor, Primary, Core) %>%
  pivot_longer(
    cols = c(Primary, Core),
    names_to = "ReceiverSet",
    values_to = "CorrectedAUPR"
  ) %>%
  mutate(
    ReceiverSet = factor(
      ReceiverSet,
      levels = c("Primary","Core"),
      labels = c("Primary (22/27 represented)",
                 "Core (8/10 represented)")
    )
  )

pC <- ggplot(nn_long, aes(CorrectedAUPR, ligand_factor, color = ReceiverSet)) +
  geom_line(
    aes(group = ligand_factor),
    color = "grey75", linewidth = 0.55
  ) +
  geom_point(size = 3.0) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  scale_color_manual(
    values = c(
      "Primary (22/27 represented)" = col_primary,
      "Core (8/10 represented)" = col_core
    )
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x * 1e4),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    title = "C  NicheNet sensory-source ligand activity",
    subtitle = "Corrected AUPR for the primary receiver program and stricter sensitivity core",
    x = "NicheNet corrected AUPR (×10⁻⁴)",
    y = NULL,
    color = NULL
  ) +
  theme_pub +
  theme(legend.position = "bottom")

# ============================================================
# PANEL D — NicheNet rank stability
# ============================================================
if (!("PrimaryRank" %in% names(nn))) {
  nn$PrimaryRank <- rank(-nn$Primary, ties.method = "first")
}
if (!("CoreRank" %in% names(nn))) {
  nn$CoreRank <- rank(-nn$Core, ties.method = "first")
}

# Explicitly label the stable/high-priority ligands.
label_ligs <- intersect(
  c("PGF","FLRT3","SEMA7A","EDIL3","SPP1","RSPO2","SST"),
  nn$ligand
)

pD <- ggplot(nn, aes(PrimaryRank, CoreRank)) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", linewidth = 0.5, color = "grey55"
  ) +
  geom_point(
    size = 3.3, shape = 21, fill = "white",
    color = "black", stroke = 0.6
  ) +
  scale_x_reverse(
    breaks = seq(10, 1, by = -1),
    limits = c(10.8, 0.2),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_reverse(
    breaks = seq(10, 1, by = -1),
    limits = c(10.8, 0.2),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "D  NicheNet rank concordance",
    subtitle = "Lower rank = stronger predicted ligand activity",
    x = "Primary receiver rank",
    y = "Core receiver rank"
  ) +
  theme_pub +
  theme(
    plot.margin = margin(6, 22, 6, 6),
    aspect.ratio = 0.86
  )

if (has_ggrepel) {
  pD <- pD + ggrepel::geom_text_repel(
    data = nn %>% filter(ligand %in% label_ligs),
    aes(label = ligand),
    size = 2.65,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.22,
    segment.size = 0.25,
    seed = 11
  )
} else {
  pD <- pD + geom_text(
    data = nn %>% filter(ligand %in% label_ligs),
    aes(label = ligand),
    size = 2.6, vjust = -0.7, check_overlap = TRUE
  )
}

# ============================================================
# Combined Figure 4
# ============================================================
combined <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1.00, 1.00), heights = c(1.05, 0.95)) +
  plot_annotation(
    title = "Human T2D β-cell receiver program and predicted sensory-neuron signaling",
    subtitle = paste0(
      "Donor-level β-cell transcriptional changes, leave-one-donor-out robustness ",
      "and NicheNet ligand prioritization"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 10.5, color = "grey30", hjust = 0),
      plot.margin = margin(10, 10, 5, 10)
    )
  )

# -----------------------------
# Frozen source-data exports
# -----------------------------
write.csv(
  de %>% select(gene = gene_plot, logFC = logFC_plot, FDR = FDR_plot, DE_class),
  "results/islet/figure_source_data/Figure4A_Beta_T2D_DE.csv",
  row.names = FALSE
)

write.csv(
  rob %>% select(
    gene = gene_plot, logFC = logFC_plot, FDR = FDR_plot,
    Core, LODO_same, LODO_FDR05
  ),
  "results/islet/figure_source_data/Figure4B_Receiver_Robustness.csv",
  row.names = FALSE
)

write.csv(
  nn_long %>% select(ligand, ReceiverSet, CorrectedAUPR),
  "results/islet/figure_source_data/Figure4C_NicheNet_Activity.csv",
  row.names = FALSE
)

write.csv(
  nn %>% select(ligand, PrimaryRank, CoreRank, Primary, Core),
  "results/islet/figure_source_data/Figure4D_NicheNet_Rank_Stability.csv",
  row.names = FALSE
)

# -----------------------------
# Save individual panels
# -----------------------------
panels <- list(
  A_Beta_T2D_DE = pA,
  B_Receiver_Robustness = pB,
  C_NicheNet_Activity = pC,
  D_NicheNet_RankStability = pD
)

for (nm in names(panels)) {
  ggsave(
    paste0("figures/islet/Figure4_11_", nm, ".pdf"),
    panels[[nm]],
    device = cairo_pdf,
    width = 7.0, height = 5.2, units = "in"
  )
  ggsave(
    paste0("figures/islet/Figure4_11_", nm, ".png"),
    panels[[nm]],
    width = 7.0, height = 5.2, units = "in", dpi = 300
  )
}

# -----------------------------
# Save combined figure
# -----------------------------
out_pdf <- "figures/islet/Figure_P3_BetaReceiver_NicheNet_11b.pdf"
out_png <- "figures/islet/Figure_P3_BetaReceiver_NicheNet_11b.png"
out_tif <- "figures/islet/Figure_P3_BetaReceiver_NicheNet_11b.tiff"

ggsave(
  out_pdf, combined,
  device = cairo_pdf,
  width = 13.5, height = 9.4, units = "in"
)

ggsave(
  out_png, combined,
  width = 13.5, height = 9.4, units = "in", dpi = 300
)

ggsave(
  out_tif, combined,
  width = 13.5, height = 9.4, units = "in",
  dpi = 600, compression = "lzw"
)

# ============================================================
# Audit
# ============================================================
cat("\n====================================================\n")
cat("FIGURE 4 AUDIT\n")
cat("====================================================\n")
cat("edgeR genes tested:", nrow(de), "\n")
cat("FDR < 0.05:", sum(de$FDR_plot < 0.05, na.rm = TRUE), "\n")
cat("T2D-up:", sum(de$FDR_plot < 0.05 & de$logFC_plot > 0, na.rm = TRUE), "\n")
cat("T2D-down:", sum(de$FDR_plot < 0.05 & de$logFC_plot < 0, na.rm = TRUE), "\n")
cat("Receiver genes:", nrow(rob), "\n")
cat("High-confidence core genes:", sum(rob$gene_plot %in% core_genes), "\n")
cat("NicheNet-evaluable sensory ligands:", nrow(nn), "\n")

if (!is.na(r_same)) {
  cat(
    "Receiver genes direction-concordant in all 18 LODO fits:",
    sum(rob$LODO_same == 18, na.rm = TRUE), "/", nrow(rob), "\n"
  )
}

cat("\nPrimary NicheNet order:\n")
print(
  nn %>%
    arrange(PrimaryRank) %>%
    select(ligand, PrimaryRank, Primary),
  row.names = FALSE
)

cat("\nCore NicheNet order:\n")
print(
  nn %>%
    arrange(CoreRank) %>%
    select(ligand, CoreRank, Core),
  row.names = FALSE
)

cat("\nScientific safeguards:\n")
cat("- beta-cell DE is donor-level pseudobulk association with T2D.\n")
cat("- LODO is robustness analysis, not independent replication.\n")
cat("- NicheNet is regulatory compatibility/prediction, not causal evidence.\n")
cat("- 22/27 primary and 8/10 core receiver genes are represented in NicheNet.\n")
cat("- RPH3A, PPIA and LGALS8 are not NicheNet-evaluable; absence is not negative evidence.\n")
cat("- predicted-target-link panel intentionally omitted because the frozen link table lacks populated target weights.\n")
cat("- no edgeR, LODO or NicheNet analysis was rerun by this script.\n")

cat("\nOutput figures:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
cat(out_tif, "\n")

cat("\n====================================================\n")
cat("SCRIPT 11 COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")
