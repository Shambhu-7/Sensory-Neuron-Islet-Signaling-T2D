#!/usr/bin/env Rscript

cat("============================================================\n")
cat("SCRIPT 18 — NATURE COMMUNICATIONS SOURCE-DATA AUDIT / STAGING\n")
cat("Paper 3: sensory neuron–islet cross-tissue manuscript\n")
cat("Inventory/staging only; NO statistical recalculation\n")
cat("============================================================\n\n")

options(stringsAsFactors = FALSE)

out_dir <- file.path("results", "source_data_package")
stage_dir <- file.path(out_dir, "staged_csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("reproducibility", recursive = TRUE, showWarnings = FALSE)

md5_one <- function(x) {
  if (is.na(x) || !file.exists(x)) return(NA_character_)
  unname(tools::md5sum(x))
}

csv_dim <- function(x) {
  if (is.na(x) || !file.exists(x) || tolower(tools::file_ext(x)) != "csv")
    return(c(NA_integer_, NA_integer_))
  z <- tryCatch(read.csv(x, check.names = FALSE), error = function(e) NULL)
  if (is.null(z)) return(c(NA_integer_, NA_integer_))
  c(nrow(z), ncol(z))
}

first_existing <- function(x) {
  y <- x[file.exists(x)]
  if (!length(y)) NA_character_ else y[1]
}

# Figure 1 is a schematic and has no numerical source-data table.
required_panels <- data.frame(
  Figure = c("Figure 1", rep("Figure 2",3), rep("Figure 3",2), "Figure 4",
             rep("Figure 5",2), rep("Supplementary Figure S1",2),
             rep("Supplementary Figure S2",3)),
  Panel = c(
    "Workflow schematic",
    "DRG differential expression", "DRG pathway remodeling", "DRG robustness",
    "Ligand-receptor candidate filtering", "Beta-cell T2D receiver",
    "NicheNet candidate activity",
    "Integrated candidate evidence", "Trace-n-Seq / subtype evidence",
    "Independent receiver validation", "Independent receptor / receiver-score validation",
    "Independent NicheNet replication", "Broad 42-ligand ranking",
    "Expression-matched empirical null"
  ),
  Numerical_Source_Data = c(FALSE, rep(TRUE,13))
)
write.csv(required_panels, file.path(out_dir,"SourceData_Required_Panels.csv"), row.names=FALSE)

# Known frozen outputs. Historical filename variants are included deliberately.
candidates <- list(
  Fig2_DRG_DE = c("results/redox/Table_P3_DE_Nageotte_vs_NearbyNeuron_edgeR.csv"),
  Fig2_DRG_Hallmark = c("results/redox/Table_P3_CAMERA_ALL_Hallmark_estimatedCorrelation.csv",
                        "results/redox/Table_P3_CAMERA_ALL_Hallmark_fixedCor001.csv"),
  Fig2_DRG_Focused = c("results/redox/Table_P3_CAMERA_Focused_Redox_estimatedCorrelation.csv",
                       "results/redox/Table_P3_CAMERA_Focused_Redox_fixedCor001.csv"),
  Fig2_DRG_LODO = c("results/redox/Table_P3_CAMERA_Focused_Robustness.csv"),
  Fig2_DRG_Anchors = c("results/redox/Table_P3_Redox_Stress_Anchor_Genes.csv"),

  Fig3_LR = c("results/islet/Table_P3_DRG_Islet_LR_Compatibility_1492cells.csv",
              "results/islet/Table_P3_DRG_Ligand_Human_Islet_Compatibility_1492cells.csv"),
  Fig3_Final13 = c("results/islet/Table_P3_TierA_SensorySource_BetaCell_HC_Shortlist.csv"),
  Fig3_Receiver = c("results/islet/Table_P3_Beta_T2D_NicheNet_PrimaryReceiverGenes.csv"),

  Fig4_Primary = c("results/nichenet/Table_P3_NicheNet_LigandActivity_Primary.csv",
                   "results/nichenet/Table_P3_NicheNet_Primary_LigandActivity.csv"),
  Fig4_Core = c("results/nichenet/Table_P3_NicheNet_LigandActivity_Core.csv",
                "results/nichenet/Table_P3_NicheNet_Core_LigandActivity.csv"),

  Fig5_A = c("results/integration/figure_source_data/Figure_Integrated_08b_PanelA_DRG.csv",
             "results/integration/figure_source_data/Figure_Integrated_PanelA_DRG.csv"),
  Fig5_B = c("results/integration/figure_source_data/Figure_Integrated_08b_PanelB_Receptors.csv",
             "results/integration/figure_source_data/Figure_Integrated_PanelB_Receptors.csv"),
  Fig5_C = c("results/integration/figure_source_data/Figure_Integrated_08b_PanelC_NicheNet.csv",
             "results/integration/figure_source_data/Figure_Integrated_PanelC_NicheNet.csv"),
  Fig5_D = c("results/integration/figure_source_data/Figure_Integrated_08b_PanelD_TraceNSeq.csv",
             "results/integration/figure_source_data/Figure_Integrated_PanelD_TraceNSeq.csv"),
  Fig5_E = c("results/integration/figure_source_data/Figure_Integrated_08b_PanelE_Subtypes.csv",
             "results/integration/figure_source_data/Figure_Integrated_PanelE_Subtypes.csv"),
  Fig5_Hierarchy = c("results/integration/Table_P3_Final_MainCandidate_Hierarchy.csv"),
  Fig5_Subtypes = c("results/integration/Table_P3_Final_SensorySubtypeEvidence.csv"),

  S1_Receiver = c("results/validation/GSE221156/Table_P3_GSE221156_Receiver_Validation.csv",
                  "results/validation/GSE221156/Table_P3_GSE221156_ReceiverValidation.csv"),
  S1_Receptors = c("results/validation/GSE221156/Table_P3_GSE221156_Receptor_Validation.csv",
                   "results/validation/GSE221156/Table_P3_GSE221156_ReceptorValidation.csv"),
  S1_Score = c("results/validation/GSE221156/Table_P3_GSE221156_ReceiverScore_ByDonor.csv",
               "results/validation/GSE221156/Table_P3_GSE221156_Receiver_Score_ByDonor.csv"),

  S2_Independent_NicheNet = c(
    "results/validation/GSE221156/NicheNet_replication/Table_P3_GSE221156_NicheNet_LigandActivity.csv",
    "results/validation/GSE221156/NicheNet_replication/Table_P3_GSE221156_Candidate_NicheNet_LigandActivity.csv"),
  S2_Broad_Discovery = c("results/nichenet/Table_P3_NicheNet_BroadHC_BetaSupported_LigandActivity.csv"),
  S2_Broad_Independent = c("results/validation/GSE221156/NicheNet_broad_replication/Table_P3_GSE221156_BroadHC_NicheNet_LigandActivity.csv")
)

# Add all final permutation CSVs as a separate audited group, because Script 16c
# already froze these outputs and exact filenames may differ by revision.
perm_dir <- "results/validation/GSE221156/NicheNet_permutation/final"
perm_csv <- if (dir.exists(perm_dir)) list.files(perm_dir, pattern="\\.csv$", full.names=TRUE) else character()
if (length(perm_csv)) candidates$S2_EmpiricalNull_FinalCSVs <- perm_csv

inventory <- do.call(rbind, lapply(names(candidates), function(k) {
  p <- candidates[[k]]
  d <- t(vapply(p, csv_dim, integer(2)))
  data.frame(Source_Key=k, Candidate_Path=p, Exists=file.exists(p),
             Size_Bytes=ifelse(file.exists(p), file.info(p)$size, NA_real_),
             MD5=vapply(p, md5_one, character(1)),
             N_Rows=d[,1], N_Columns=d[,2], check.names=FALSE)
}))
write.csv(inventory, file.path(out_dir,"SourceData_File_Inventory.csv"), row.names=FALSE)

selected <- data.frame(
  Source_Key=names(candidates),
  Selected_Path=vapply(candidates, first_existing, character(1)),
  stringsAsFactors=FALSE
)
selected$Found <- !is.na(selected$Selected_Path)
selected$MD5 <- vapply(selected$Selected_Path, md5_one, character(1))
dims <- t(vapply(selected$Selected_Path, csv_dim, integer(2)))
selected$N_Rows <- dims[,1]
selected$N_Columns <- dims[,2]

selected$Staged_File <- NA_character_
for (i in seq_len(nrow(selected))) {
  src <- selected$Selected_Path[i]
  if (!is.na(src) && file.exists(src) && tolower(tools::file_ext(src)) == "csv") {
    nm <- sprintf("%02d_%s.csv", i, gsub("[^A-Za-z0-9_]+","_",selected$Source_Key[i]))
    dst <- file.path(stage_dir,nm)
    if (file.copy(src,dst,overwrite=TRUE)) selected$Staged_File[i] <- dst
  }
}
write.csv(selected, file.path(out_dir,"SourceData_Staged_Files.csv"), row.names=FALSE)

# Broad discovery of likely frozen source tables helps recover unmapped filenames.
all_csv <- if (dir.exists("results")) list.files("results", pattern="\\.csv$", recursive=TRUE, full.names=TRUE) else character()
likely <- all_csv[grepl("P3|GSE295206|GSE81608|GSE221156|NicheNet|Trace|Candidate|Receiver|Ligand|Receptor|Integrated|Source", all_csv, ignore.case=TRUE)]
discovered <- data.frame(Path=likely,
                         Size_Bytes=if(length(likely)) file.info(likely)$size else numeric(),
                         MD5=if(length(likely)) vapply(likely,md5_one,character(1)) else character())
write.csv(discovered, file.path(out_dir,"SourceData_Discovered_Likely_CSVs.csv"), row.names=FALSE)

# Panel-level readiness. These checks concern file mapping only, not scientific validity.
has <- setNames(selected$Found, selected$Source_Key)
ready <- function(keys) all(has[keys] %in% TRUE)
panel_audit <- data.frame(
  Figure=c("Figure 2","Figure 2","Figure 2","Figure 3","Figure 3","Figure 4","Figure 5","Figure 5","Supplementary Figure S1","Supplementary Figure S1","Supplementary Figure S2","Supplementary Figure S2","Supplementary Figure S2"),
  Component=c("DRG differential expression","DRG pathway remodeling","DRG robustness","Ligand-receptor candidate filtering","Beta-cell receiver","NicheNet candidate activity","Integrated five-candidate evidence","Trace-n-Seq / subtype evidence","Independent receiver validation","Independent receptor / receiver-score validation","Independent NicheNet replication","Broad 42-ligand cross-cohort ranking","Expression-matched empirical null"),
  Ready=c(
    ready("Fig2_DRG_DE"),
    ready(c("Fig2_DRG_Hallmark","Fig2_DRG_Focused")),
    ready("Fig2_DRG_LODO"),
    ready(c("Fig3_LR","Fig3_Final13")),
    ready("Fig3_Receiver"),
    ready(c("Fig4_Primary","Fig4_Core")),
    ready(c("Fig5_A","Fig5_B","Fig5_C","Fig5_Hierarchy")),
    ready(c("Fig5_D","Fig5_E","Fig5_Subtypes")),
    ready("S1_Receiver"),
    ready(c("S1_Receptors","S1_Score")),
    ready("S2_Independent_NicheNet"),
    ready(c("S2_Broad_Discovery","S2_Broad_Independent")),
    "S2_EmpiricalNull_FinalCSVs" %in% names(has) && isTRUE(has["S2_EmpiricalNull_FinalCSVs"])
  ), stringsAsFactors=FALSE
)
panel_audit$Status <- ifelse(panel_audit$Ready,"READY_FOR_ASSEMBLY","NEEDS_FILENAME_MAPPING")
write.csv(panel_audit, file.path(out_dir,"SourceData_Panel_Candidate_Audit.csv"), row.names=FALSE)
write.csv(panel_audit[!panel_audit$Ready,,drop=FALSE],
          file.path(out_dir,"SourceData_Panels_Needing_Mapping.csv"), row.names=FALSE)

cat("\n============================================================\n")
cat("SOURCE-DATA AUDIT\n")
cat("============================================================\n")
cat("\nFigure 1: workflow schematic — no numerical source-data table required.\n")
cat("\nSELECTED FROZEN FILES\n")
print(selected[,c("Source_Key","Found","Selected_Path","N_Rows","N_Columns")], row.names=FALSE)
cat("\nPANEL READINESS\n")
print(panel_audit[,c("Figure","Component","Status")], row.names=FALSE)
cat("\nSUMMARY\n")
cat("Source groups audited:", nrow(selected), "\n")
cat("Source groups found:", sum(selected$Found), "\n")
cat("Source groups needing filename mapping:", sum(!selected$Found), "\n")
cat("Panels ready:", sum(panel_audit$Ready), "/", nrow(panel_audit), "\n")
cat("Likely CSVs discovered:", nrow(discovered), "\n")
cat("\nOutputs:", out_dir, "\n")
cat("No statistics were recalculated.\n")

writeLines(capture.output(sessionInfo()), "reproducibility/18_SourceData_audit_sessionInfo.txt")

cat("\n============================================================\n")
cat("SCRIPT 18 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
