cat("====================================================\n")
cat("SCRIPT 08: INTEGRATED CANDIDATE FIGURE\n")
cat("====================================================\n\n")

required <- c("ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse=", "))
library(ggplot2)
library(patchwork)

dir.create("figures/integration", recursive=TRUE, showWarnings=FALSE)
dir.create("results/integration/figure_source_data", recursive=TRUE, showWarnings=FALSE)
dir.create("reproducibility", recursive=TRUE, showWarnings=FALSE)

f_main <- "results/integration/Table_P3_Final_MainCandidate_Hierarchy.csv"
f_sub <- "results/integration/Table_P3_Final_SensorySubtypeEvidence.csv"
if (!file.exists(f_main)) stop("Missing: ", f_main)
if (!file.exists(f_sub)) stop("Missing: ", f_sub)

x <- read.csv(f_main, check.names=FALSE, stringsAsFactors=FALSE)
sub <- read.csv(f_sub, check.names=FALSE, stringsAsFactors=FALSE)
genes <- c("SEMA7A","EDIL3","FLRT3","PGF","SPP1")
if (!all(genes %in% x$HumanGene)) stop("Frozen five-candidate hierarchy is incomplete.")
x <- x[match(genes, x$HumanGene), , drop=FALSE]
lev <- rev(genes)
x$HumanGene <- factor(x$HumanGene, levels=lev)

theme_paper <- theme_classic(base_size=11) + theme(
  axis.title=element_text(face="bold"), axis.text=element_text(color="black"),
  plot.title=element_text(face="bold", size=12, hjust=0),
  plot.subtitle=element_text(size=9), legend.title=element_text(face="bold")
)

# A: human diabetic DRG
A <- data.frame(HumanGene=x$HumanGene,
                NearbyNeuron_Enrichment=-x$LR_logFC,
                FDR=x$LR_FDR, Tier=x$LR_candidate_tier)
A$MinusLog10FDR <- -log10(pmax(A$FDR, .Machine$double.xmin))
write.csv(A,"results/integration/figure_source_data/Figure_Integrated_PanelA_DRG.csv",row.names=FALSE)
pA <- ggplot(A,aes(NearbyNeuron_Enrichment,HumanGene)) +
  geom_segment(aes(x=0,xend=NearbyNeuron_Enrichment,yend=HumanGene),linewidth=.7) +
  geom_point(aes(size=MinusLog10FDR,shape=Tier)) +
  scale_size_continuous(name=expression(-log[10]~FDR),range=c(3,7)) +
  labs(title="A  Human diabetic DRG",subtitle="Nearby-neuron-associated ligand enrichment",
       x=expression("Nearby-neuron enrichment (" * -log[2] * "FC)"),y=NULL,shape="DRG tier") +
  coord_cartesian(xlim=c(0,max(A$NearbyNeuron_Enrichment,na.rm=TRUE)*1.12)) + theme_paper

# B: beta-cell receptor compatibility
B <- rbind(
  data.frame(HumanGene=as.character(x$HumanGene),Support="Supported in ND + T2D",
             N_receptors=x$LR_n_HC_beta_both_receptors),
  data.frame(HumanGene=as.character(x$HumanGene),Support="Strong in ND + T2D",
             N_receptors=x$LR_n_HC_beta_strong_both_receptors)
)
B$HumanGene <- factor(B$HumanGene,levels=lev)
B$Support <- factor(B$Support,levels=c("Supported in ND + T2D","Strong in ND + T2D"))
write.csv(B,"results/integration/figure_source_data/Figure_Integrated_PanelB_Receptors.csv",row.names=FALSE)
pB <- ggplot(B,aes(N_receptors,HumanGene,fill=Support)) +
  geom_col(position=position_dodge(width=.72),width=.62) +
  scale_x_continuous(breaks=0:5,expand=expansion(mult=c(0,.08))) +
  labs(title="B  Human beta-cell receptors",subtitle="High-confidence receptor compatibility",
       x="Number of beta-cell receptors",y=NULL,fill=NULL) + theme_paper +
  theme(legend.position="bottom")

# C: NicheNet
C <- rbind(
  data.frame(HumanGene=as.character(x$HumanGene),Receiver="Primary 22/27 receiver genes",
             Corrected_AUPR=x$NicheNet_Primary_AUPR_corrected,Rank=x$NicheNet_Primary_Rank),
  data.frame(HumanGene=as.character(x$HumanGene),Receiver="Core 8/10 receiver genes",
             Corrected_AUPR=x$NicheNet_Core_AUPR_corrected,Rank=x$NicheNet_Core_Rank)
)
C$HumanGene <- factor(C$HumanGene,levels=lev)
C$Receiver <- factor(C$Receiver,levels=c("Primary 22/27 receiver genes","Core 8/10 receiver genes"))
write.csv(C,"results/integration/figure_source_data/Figure_Integrated_PanelC_NicheNet.csv",row.names=FALSE)
pC <- ggplot(C,aes(Corrected_AUPR,HumanGene)) +
  geom_vline(xintercept=0,linetype="dashed",linewidth=.4) +
  geom_line(aes(group=HumanGene),linewidth=.6) +
  geom_point(aes(shape=Receiver),size=3) +
  geom_text(aes(label=paste0("#",Rank)),nudge_x=max(C$Corrected_AUPR,na.rm=TRUE)*.035,
            hjust=0,size=3,show.legend=FALSE) +
  scale_x_continuous(expand=expansion(mult=c(.04,.18))) +
  labs(title="C  Beta-cell receiver prediction",subtitle="NicheNet ligand activity",
       x="Corrected AUPR",y=NULL,shape=NULL) + theme_paper + theme(legend.position="bottom")

# D: pancreas-traced DRG and NG
D <- rbind(
 data.frame(HumanGene=as.character(x$HumanGene),Cohort="DRG batch 1",PctDetected=x$DRG_Batch1_PctDetected),
 data.frame(HumanGene=as.character(x$HumanGene),Cohort="DRG batch 2",PctDetected=x$DRG_Batch2_PctDetected),
 data.frame(HumanGene=as.character(x$HumanGene),Cohort="NG",PctDetected=x$NG_PctDetected)
)
D$HumanGene <- factor(D$HumanGene,levels=lev)
D$Cohort <- factor(D$Cohort,levels=c("DRG batch 1","DRG batch 2","NG"))
write.csv(D,"results/integration/figure_source_data/Figure_Integrated_PanelD_TraceNSeq.csv",row.names=FALSE)
pD <- ggplot(D,aes(Cohort,HumanGene,fill=PctDetected)) + geom_tile(linewidth=.5) +
  geom_text(aes(label=sprintf("%.0f",PctDetected)),size=3.3) +
  scale_fill_viridis_c(limits=c(0,100),name="% detected") +
  labs(title="D  Pancreas-traced neurons",subtitle="Trace-n-Seq anatomical validation",x=NULL,y=NULL) +
  theme_paper + theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="bottom")

# E: sensory subtype evidence
need <- c("HumanGene","CellType","N_cells","Pct_detected","Mean_log1p_CPM")
if (!all(need %in% names(sub))) stop("Subtype table lacks expected columns: ",paste(setdiff(need,names(sub)),collapse=", "))
sub <- sub[sub$HumanGene %in% genes,,drop=FALSE]
short <- c("Peripheral sensory neurofilament neurons"="Neurofilament",
           "Peripheral sensory non-peptidergic neurons"="Non-peptidergic",
           "Peripheral sensory peptidergic neurons"="Peptidergic")
sub$SensorySubtype <- unname(short[as.character(sub$CellType)])
sub <- sub[!is.na(sub$SensorySubtype),,drop=FALSE]
sub$HumanGene <- factor(sub$HumanGene,levels=lev)
sub$SensorySubtype <- factor(sub$SensorySubtype,levels=c("Neurofilament","Non-peptidergic","Peptidergic"))
write.csv(sub,"results/integration/figure_source_data/Figure_Integrated_PanelE_Subtypes.csv",row.names=FALSE)
pE <- ggplot(sub,aes(SensorySubtype,HumanGene,fill=Pct_detected)) + geom_tile(linewidth=.5) +
  geom_text(aes(label=sprintf("%.0f",Pct_detected)),size=3.3) +
  scale_fill_viridis_c(limits=c(0,100),name="% detected") +
  labs(title="E  Sensory-neuron subtype",subtitle="Combined pancreas-traced DRG cohorts",x=NULL,y=NULL) +
  theme_paper + theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="bottom")

combined <- ((pA+pB+pC)+plot_layout(widths=c(1.15,1,1.15))) /
            ((pD+pE)+plot_layout(widths=c(1,1))) +
  plot_annotation(title="Integrated prioritization of candidate sensory neuron-islet signaling axes",
                  subtitle="Human diabetic DRG spatial evidence, human beta-cell receptor compatibility, NicheNet receiver prediction and pancreas-traced sensory-neuron expression",
                  theme=theme(plot.title=element_text(face="bold",size=15),plot.subtitle=element_text(size=10)))

panels <- list(A_DRG=pA,B_Receptors=pB,C_NicheNet=pC,D_TraceNSeq=pD,E_Subtypes=pE)
for (nm in names(panels)) {
 ggsave(paste0("figures/integration/Figure_Integrated_",nm,".pdf"),panels[[nm]],width=6.2,height=4.2,units="in")
 ggsave(paste0("figures/integration/Figure_Integrated_",nm,".png"),panels[[nm]],width=6.2,height=4.2,units="in",dpi=400)
}
ggsave("figures/integration/Figure_P3_IntegratedCandidateEvidence.pdf",combined,width=13,height=9,units="in")
ggsave("figures/integration/Figure_P3_IntegratedCandidateEvidence.png",combined,width=13,height=9,units="in",dpi=400)
ggsave("figures/integration/Figure_P3_IntegratedCandidateEvidence.tiff",combined,width=13,height=9,units="in",dpi=600,compression="lzw")

cat("\n====================================================\nFIGURE SOURCE-DATA AUDIT\n====================================================\n")
cat("Panel A candidates:",nrow(A),"\nPanel B rows:",nrow(B),"\nPanel C rows:",nrow(C),"\nPanel D rows:",nrow(D),"\nPanel E rows:",nrow(sub),"\n")
cat("Candidate order:",paste(genes,collapse=" -> "),"\n")
cat("\nIMPORTANT:\nNegative original DRG logFC denotes nearby-neuron enrichment.\nPanel A displays -logFC as positive nearby-neuron enrichment.\nReceptor, Trace-n-Seq and subtype panels are descriptive.\nNicheNet ranks are model-derived ligand-activity ranks.\nNo composite statistical score is calculated.\n")
writeLines(capture.output(sessionInfo()),"reproducibility/08_integrated_candidate_figure_sessionInfo.txt")
cat("\nOutput figures:\n  figures/integration/Figure_P3_IntegratedCandidateEvidence.pdf\n  figures/integration/Figure_P3_IntegratedCandidateEvidence.tiff\n  figures/integration/Figure_P3_IntegratedCandidateEvidence.png\n")
cat("\n====================================================\nSCRIPT 08 COMPLETED SUCCESSFULLY\n====================================================\n")
