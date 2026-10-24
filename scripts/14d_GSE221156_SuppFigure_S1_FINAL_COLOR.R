#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

cat(strrep("=", 78), "\n", sep="")
cat("SUPPLEMENTARY FIGURE S1 — INDEPENDENT GSE221156 BETA-CELL VALIDATION\n")
cat(strrep("=", 78), "\n", sep="")

indir <- "results/validation/GSE221156"
outdir <- "figures/supplementary"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

f_meta <- file.path(indir, "GSE221156_Donor_Metadata_Audit.csv")
f_val  <- file.path(indir, "Table_P3_GSE221156_27Gene_Receiver_Validation.csv")
f_sum  <- file.path(indir, "Table_P3_GSE221156_ReceiverProgram_ValidationSummary.csv")
f_cam  <- file.path(indir, "Table_P3_GSE221156_ReceiverProgram_CAMERA_PR.csv")
f_down <- file.path(indir, "Table_P3_GSE221156_ReceiverProgram_DownGeneAudit.csv")
f_lr   <- file.path(indir, "Table_P3_GSE221156_Ligand_Receptor_Reproduction.csv")
f_score<- file.path(indir, "Table_P3_GSE221156_ReceiverScore_All48Donors.csv")

req <- c(f_meta,f_val,f_sum,f_cam,f_down,f_lr,f_score)
miss <- req[!file.exists(req)]
if(length(miss)) stop("Missing required file(s):\n", paste(miss, collapse="\n"))
cat("INPUT FILES: PASS\n")

meta <- read.csv(f_meta, check.names=FALSE)
val <- read.csv(f_val, check.names=FALSE)
sumtab <- read.csv(f_sum, check.names=FALSE)
cam <- read.csv(f_cam, check.names=FALSE)
down <- read.csv(f_down, check.names=FALSE)
lr <- read.csv(f_lr, check.names=FALSE)
score <- read.csv(f_score, check.names=FALSE)

# -----------------------------------------------------------------------------
# Global publication theme
# -----------------------------------------------------------------------------
# Restrained semantic palette used consistently across panels
disease_cols <- c(
  "Non-diabetic" = "#4C78A8",
  "Prediabetes"  = "#E3A33A",
  "T2D"          = "#D95F59"
)
replication_cols <- c(
  "Same direction" = "#2A9D8F",
  "Opposite direction" = "#C44E52"
)
test_cols <- c(
  "Directional concordance" = "#2A9D8F",
  "CAMERA-PR enrichment" = "#7B4AB5"
)

theme_pub <- theme_classic(base_size=10) +
  theme(
    plot.title = element_text(face="bold", size=11, margin=margin(b=4)),
    plot.subtitle = element_text(size=8.7, colour="grey25", margin=margin(b=7)),
    axis.title = element_text(size=9.2),
    axis.text = element_text(size=8.5, colour="black"),
    legend.title = element_text(size=8.5, face="bold"),
    legend.text = element_text(size=8.2),
    plot.margin = margin(7,8,7,8)
  )

# -----------------------------------------------------------------------------
# A. Independent cohort design
# -----------------------------------------------------------------------------
disease_levels <- c("normal","prediabetes syndrome","type 2 diabetes mellitus")
disease_labels <- c("Non-diabetic","Prediabetes","T2D")
meta$disease <- factor(meta$disease, levels=disease_levels)
donor_counts <- table(meta$disease)
cell_counts <- tapply(meta$n_beta_cells, meta$disease, sum)
cohort_df <- data.frame(
  disease=factor(disease_levels, levels=disease_levels),
  donors=as.integer(donor_counts[disease_levels]),
  cells=as.numeric(cell_counts[disease_levels])
)
cohort_df$label <- factor(disease_labels, levels=disease_labels)

pA <- ggplot(cohort_df, aes(x=label, y=donors, fill=label)) +
  geom_col(width=.62, alpha=.88, show.legend=FALSE) +
  scale_fill_manual(values=disease_cols) +
  geom_text(aes(label=donors), vjust=-.35, size=3.3, fontface="bold") +
  annotate("text", x=2, y=max(cohort_df$donors)+2.4,
           label=paste0(format(sum(cohort_df$cells), big.mark=","), " β cells"),
           size=3.2, fontface="bold") +
  scale_y_continuous(expand=expansion(mult=c(0,.18)), breaks=seq(0,20,5)) +
  labs(title="Independent human β-cell validation cohort",
       subtitle="GSE221156 · 48 donors · donor-level analysis",
       x=NULL, y="Donors") + theme_pub +
  theme(axis.text.x=element_text(face="bold"))

# -----------------------------------------------------------------------------
# B. Original vs independent validation effect sizes
# -----------------------------------------------------------------------------
v <- val[val$tested_by_edgeR %in% TRUE & is.finite(val$primary_logFC) & is.finite(val$logFC),,drop=FALSE]
if(nrow(v) < 3) stop("Too few receiver genes available for panel B.")
rho <- suppressWarnings(cor(v$primary_logFC, v$logFC, method="spearman"))
rhop <- suppressWarnings(cor.test(v$primary_logFC, v$logFC, method="spearman", exact=FALSE)$p.value)
k <- sum(v$same_direction %in% TRUE); n <- nrow(v)
bp <- binom.test(k,n,p=.5,alternative="greater")$p.value
v$replication <- ifelse(v$same_direction %in% TRUE, "Same direction", "Opposite direction")
v$label <- ifelse(v$gene %in% c("CLU","TTR","NPNT"), v$gene, "")
lim <- max(abs(c(v$primary_logFC,v$logFC)), na.rm=TRUE) * 1.08

pB <- ggplot(v, aes(primary_logFC, logFC)) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey75") +
  geom_vline(xintercept=0, linewidth=.3, colour="grey75") +
  geom_abline(slope=1, intercept=0, linetype=2, linewidth=.45, colour="grey55") +
  geom_point(aes(shape=replication, colour=replication), size=2.8, stroke=.85, fill="white") +
  geom_text(aes(label=label, colour=replication), check_overlap=TRUE, nudge_y=.12, size=2.7, fontface="bold", show.legend=FALSE) +
  scale_shape_manual(values=c("Same direction"=21,"Opposite direction"=4)) +
  scale_colour_manual(values=replication_cols) +
  coord_cartesian(xlim=c(-lim,lim), ylim=c(-lim,lim), clip="off") +
  labs(title="Receiver-program replication",
       subtitle=sprintf("%d/%d same direction · Spearman ρ = %.2f, P = %.2g", k,n,rho,rhop),
       x="GSE81608 T2D log2FC", y="GSE221156 T2D log2FC", shape=NULL) +
  theme_pub + theme(legend.position="bottom")

# -----------------------------------------------------------------------------
# C. Program-level statistical validation
# -----------------------------------------------------------------------------
# Combine prespecified directional-concordance and CAMERA-PR tests in one panel.
ps <- sumtab
ps$short <- ifelse(grepl("Primary", ps$gene_set), "27-gene receiver", "10-gene core")
cam$short <- ifelse(grepl("Original_receiver", cam$gene_set), "27-gene receiver", "10-gene core")

c1 <- data.frame(
  set=ps$short,
  test="Directional concordance",
  P=ps$exact_binomial_one_sided_P,
  detail=sprintf("%d/%d", ps$n_same_direction, ps$n_tested_edgeR)
)
c2 <- data.frame(
  set=cam$short,
  test="CAMERA-PR enrichment",
  P=cam$PValue,
  detail=paste0("FDR=", formatC(cam$FDR, format="g", digits=2))
)
cdf <- rbind(c1,c2)
cdf$score <- -log10(cdf$P)
cdf$set <- factor(cdf$set, levels=c("10-gene core","27-gene receiver"))
cdf$test <- factor(cdf$test, levels=c("Directional concordance","CAMERA-PR enrichment"))

cdf$ybase <- ifelse(cdf$set == "27-gene receiver", 2, 1)
cdf$y <- cdf$ybase + ifelse(cdf$test == "Directional concordance", -0.11, 0.11)

pC <- ggplot(cdf, aes(x=score, y=y, shape=test, colour=test)) +
  geom_vline(xintercept=-log10(.05), linetype=2, linewidth=.45, colour="grey55") +
  geom_segment(aes(x=0, xend=score, yend=y, colour=test), linewidth=.7, alpha=.75) +
  geom_point(size=3.2, fill="white", stroke=.9) +
  geom_text(aes(label=detail, colour=test), hjust=-.15, vjust=.35, size=2.65, show.legend=FALSE) +
  scale_shape_manual(values=c("Directional concordance"=21,"CAMERA-PR enrichment"=22)) +
  scale_colour_manual(values=test_cols) +
  scale_y_continuous(breaks=c(1,2), labels=c("10-gene core","27-gene receiver"),
                     limits=c(.65,2.35)) +
  scale_x_continuous(expand=expansion(mult=c(0,.38))) +
  labs(title="Prespecified program-level validation",
       subtitle="Dashed line: P = 0.05; CAMERA tests original T2D-up program",
       x=expression(-log[10](P)), y=NULL, shape=NULL) + theme_pub +
  theme(legend.position="bottom")

# -----------------------------------------------------------------------------
# D. Independent receptor reproduction + disease-stage receiver score
# -----------------------------------------------------------------------------
reproduced <- sum(lr$beta_receptor_support_reproduced %in% TRUE)
total_lr <- nrow(lr)

score$disease <- factor(score$disease, levels=disease_levels, labels=disease_labels)
med <- aggregate(receiver_direction_aligned_score ~ disease, score, median)
medlab <- paste(sprintf("%s %.2f", med$disease, med$receiver_direction_aligned_score), collapse="  ·  ")

pD <- ggplot(score, aes(x=disease, y=receiver_direction_aligned_score)) +
  geom_hline(yintercept=0, linewidth=.35, colour="grey70") +
  geom_boxplot(aes(fill=disease), width=.58, outlier.shape=NA, colour="grey20", alpha=.28, linewidth=.6, show.legend=FALSE) +
  geom_jitter(aes(colour=disease), width=.12, height=0, size=1.55, alpha=.78, show.legend=FALSE) +
  scale_fill_manual(values=disease_cols) +
  scale_colour_manual(values=disease_cols) +
  labs(
    title="Receptor reproduction and direction-aligned receiver score",
    subtitle=paste0(reproduced,"/",total_lr,
                    " candidate ligand axes retain β-cell receptor support in both ND and T2D  ·  ",
                    medlab),
    x=NULL, y="Receiver score"
  ) + theme_pub +
  theme(axis.text.x=element_text(size=8.2, face="bold"),
        plot.subtitle=element_text(size=8.3),
        plot.margin=margin(7,10,7,7))

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title="Independent human β-cell validation in GSE221156",
    subtitle="Prespecified replication of the GSE81608 T2D receiver program and β-cell receptor compatibility",
    tag_levels="A",
    theme=theme(
      plot.title=element_text(face="bold",size=15),
      plot.subtitle=element_text(size=10.5,colour="grey25"),
      plot.tag=element_text(face="bold",size=13)
    )
  )

png_file <- file.path(outdir,"Figure_S1_GSE221156_Independent_BetaCell_Validation_FINAL_COLOR.png")
pdf_file <- file.path(outdir,"Figure_S1_GSE221156_Independent_BetaCell_Validation_FINAL_COLOR.pdf")

ggsave(png_file, fig, width=15.5, height=9.0, dpi=600, bg="white")
ggsave(pdf_file, fig, width=15.5, height=9.0, device=cairo_pdf, bg="white")

cat("\nKEY VALUES USED\n")
cat(sprintf("  Cohort: %d donors; %s beta cells\n", nrow(meta), format(sum(meta$n_beta_cells),big.mark=",")))
cat(sprintf("  Receiver: %d/%d same direction; rho=%.4f; P=%.6g; binomial P=%.6g\n",k,n,rho,rhop,bp))
cat(sprintf("  Receptor compatibility reproduced: %d/%d candidate ligand axes\n",reproduced,total_lr))
if(nrow(down)) {
  for(i in seq_len(nrow(down))) cat(sprintf("  %s validation logFC %.4f; FDR %.6g\n",down$gene[i],down$logFC[i],down$FDR[i]))
}
cat("\nOUTPUTS\n")
cat("  ",png_file,"\n",sep="")
cat("  ",pdf_file,"\n",sep="")
cat("\nINTERPRETATION\n")
cat("  Independent beta-cell transcriptional and receptor validation only.\n")
cat("  The normal-to-prediabetes-to-T2D receiver score is descriptive, not a formal trajectory test.\n")
cat("  No direct sensory-neuron-to-beta-cell signaling or causality is inferred.\n")
cat("\n",strrep("=",78),"\n",sep="")
cat("SUPPLEMENTARY FIGURE S1 COMPLETE\n")
cat(strrep("=",78),"\n",sep="")
