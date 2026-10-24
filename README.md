# Paper 3 reproducibility archive

## Manuscript
**Integrated human transcriptomics prioritizes candidate sensory neuron–islet signaling axes in type 2 diabetes**

Author: Shambhu Yadav, PhD  
Department of Biotechnology, Manipal Institute of Technology, Manipal Academy of Higher Education, Manipal, Karnataka, India  
ORCID: 0000-0002-0622-1854

## Scope
This repository contains the analysis scripts, frozen tabular outputs, reproducibility records and journal Source Data workbook used for the manuscript. It is a computational reproducibility archive. It does **not** establish direct sensory neuron–beta-cell signaling or causality.

## Public datasets
- **GSE295206** — human diabetic dorsal root ganglion spatial transcriptomics; discovery DRG analysis.
- **GSE81608** — human pancreatic islet single-cell RNA-seq; primary beta-cell receptor compatibility and T2D receiver analysis. The published 1,492-cell representation is primary; the current 1,600-cell GEO representation is representation sensitivity, not an independent cohort.
- **E-MTAB-12940** — mouse Trace-n-Seq pancreas-traced sensory neurons; anatomical-expression support only.
- **GSE221156** — independent human beta-cell validation cohort (48 donors; 99,029 beta cells).

Raw public datasets are intentionally not redistributed. Download them from the source repositories using the accessions above and the included download/audit scripts.

## Large-file provenance
- GSE221156 beta-cell H5AD: 2,477,777,480 bytes; MD5 `e9097864a452e6983f35dce8aafad262`.
- E-MTAB-12940 `smartseq_raw_counts.txt`: 174,315,321 bytes; MD5 `39d43906638b2ae3d87009905c49806f`.
- E-MTAB-12940 `smartseq_raw_counts2.txt`: 182,294,665 bytes; MD5 `02a141ff535d351963fd4e8baa3a8629`.
- Frozen NicheNet ligand-target prior MD5: `b09606b04b2d4490418d9028c0e58b9f`.

## Directory structure
- `scripts/` — analysis, validation and figure-generation scripts.
- `results/` — frozen tabular outputs used by the manuscript and figures.
- `reproducibility/` — console logs, session information and audits.
- `source_data/` — Nature Communications Source Data workbook, when present.
- `provenance/` — selected retrieval/model metadata.
- `MANIFEST_SHA256.csv` — SHA-256 checksum and size for every archived file.

## Key frozen validation results
Independent broad-universe NicheNet ranks were concordant between GSE81608 and GSE221156 (Spearman rho = 0.548, P = 1.71e-4). Expression-matched empirical-null sensitivity testing used 10,000 matched random receiver programs per cohort; no ligand achieved empirical FDR < 0.05 in both cohorts. These results are interpreted as sensitivity analyses and regulatory compatibility, not direct signaling evidence.

## Reproducibility notes
Scripts are numbered by workflow stage. Preserve the frozen thresholds, candidate universe, receiver definitions and NicheNet prior described in the manuscript. Do not substitute the GSE81608 1,600-cell representation for the primary 1,492-cell analysis.

## Software environment
Core analyses were performed in R 4.5.2 on arm64 macOS. Key packages include edgeR 4.8.2, limma 3.66.0, nichenetr 2.2.0, ggplot2 4.0.1 and patchwork 1.3.2. Python 3.14.3 was used for selected download/audit tasks. Exact session records are provided under `reproducibility/` where generated.

## License
Code in `scripts/` is released under the MIT License. Public-source data and derived outputs remain subject to the terms of their originating repositories/publications.
