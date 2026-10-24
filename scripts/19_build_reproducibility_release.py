#!/usr/bin/env python3
from pathlib import Path
import shutil, hashlib, csv, zipfile, os, sys, json
from datetime import date

ROOT = Path.home() / "Desktop/RNA3.0/Paper3_Sensory_Islet"
if not ROOT.exists():
    raise SystemExit(f"Project root not found: {ROOT}")

OUT = ROOT / "release_v1.0.0"
GITHUB = OUT / "github_repository"
ZENODO = OUT / "zenodo_archive"

# Fresh deterministic staging area.
if OUT.exists():
    shutil.rmtree(OUT)
GITHUB.mkdir(parents=True)
ZENODO.mkdir(parents=True)

print("="*72)
print("SCRIPT 19 — PAPER 3 REPRODUCIBILITY RELEASE BUILDER")
print("Inventory/package only; NO statistical recalculation")
print("="*72)

# ---------- helpers ----------
def sha256(path, block=1024*1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b: break
            h.update(b)
    return h.hexdigest()

def copy_file(src, dst):
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

def copy_tree_filtered(src, dst, allowed_ext=None, max_bytes=None, exclude_parts=()):
    src, dst = Path(src), Path(dst)
    if not src.exists():
        return 0
    n = 0
    for p in sorted(src.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(src)
        if any(part in exclude_parts for part in rel.parts):
            continue
        if allowed_ext and p.suffix.lower() not in allowed_ext:
            continue
        if max_bytes is not None and p.stat().st_size > max_bytes:
            continue
        copy_file(p, dst / rel)
        n += 1
    return n

def make_zip(src_dir, zip_path):
    if zip_path.exists(): zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(src_dir.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(src_dir.parent))

# ---------- Source Data workbook ----------
workbook_candidates = [
    ROOT / "Yadav_Paper3_Source_Data.xlsx",
    Path.home() / "Downloads/Yadav_Paper3_Source_Data.xlsx",
]
workbook = next((p for p in workbook_candidates if p.exists()), None)
if workbook:
    print("FOUND Source Data workbook:", workbook)
else:
    print("WARNING: Yadav_Paper3_Source_Data.xlsx not found in project root or Downloads.")
    print("         Release will still be built; copy the workbook into the project and rerun before final publication.")

# ---------- GitHub package ----------
# Code: preserve all numbered/final scripts for auditability. Raw data are excluded.
script_count = copy_tree_filtered(
    ROOT / "scripts", GITHUB / "scripts",
    allowed_ext={".r", ".R", ".py", ".sh"}, max_bytes=20*1024*1024
)

# Reproducibility records (text/log/session info only).
repro_count = copy_tree_filtered(
    ROOT / "reproducibility", GITHUB / "reproducibility",
    allowed_ext={".txt", ".log", ".csv", ".json", ".md"}, max_bytes=20*1024*1024
)

# Frozen tabular results. Exclude bulky binary matrices and raw data.
result_count = copy_tree_filtered(
    ROOT / "results", GITHUB / "results",
    allowed_ext={".csv", ".tsv", ".txt", ".json", ".md"}, max_bytes=50*1024*1024,
    exclude_parts={"smoketest"}
)

if workbook:
    copy_file(workbook, GITHUB / "source_data" / workbook.name)

# Public-data provenance/checksum records if present.
for rel in [
    "data/raw/GSE221156/GSE221156_CELLxGENE_download_record.json",
    "results/nichenet/NicheNet_Model_and_Analysis_Record.csv",
    "results/ligands/OmniPath_DIRECT_retrieval_metadata.csv",
    "results/ligands/OmniPath_FROZEN_LR_network_metadata.csv",
]:
    p = ROOT / rel
    if p.exists(): copy_file(p, GITHUB / "provenance" / Path(rel).name)

readme = '''# Paper 3 reproducibility archive\n\n## Manuscript\n**Integrated human transcriptomics prioritizes candidate sensory neuron–islet signaling axes in type 2 diabetes**\n\nAuthor: Shambhu Yadav, PhD  \nDepartment of Biotechnology, Manipal Institute of Technology, Manipal Academy of Higher Education, Manipal, Karnataka, India  \nORCID: 0000-0002-0622-1854\n\n## Scope\nThis repository contains the analysis scripts, frozen tabular outputs, reproducibility records and journal Source Data workbook used for the manuscript. It is a computational reproducibility archive. It does **not** establish direct sensory neuron–beta-cell signaling or causality.\n\n## Public datasets\n- **GSE295206** — human diabetic dorsal root ganglion spatial transcriptomics; discovery DRG analysis.\n- **GSE81608** — human pancreatic islet single-cell RNA-seq; primary beta-cell receptor compatibility and T2D receiver analysis. The published 1,492-cell representation is primary; the current 1,600-cell GEO representation is representation sensitivity, not an independent cohort.\n- **E-MTAB-12940** — mouse Trace-n-Seq pancreas-traced sensory neurons; anatomical-expression support only.\n- **GSE221156** — independent human beta-cell validation cohort (48 donors; 99,029 beta cells).\n\nRaw public datasets are intentionally not redistributed. Download them from the source repositories using the accessions above and the included download/audit scripts.\n\n## Large-file provenance\n- GSE221156 beta-cell H5AD: 2,477,777,480 bytes; MD5 `e9097864a452e6983f35dce8aafad262`.\n- E-MTAB-12940 `smartseq_raw_counts.txt`: 174,315,321 bytes; MD5 `39d43906638b2ae3d87009905c49806f`.\n- E-MTAB-12940 `smartseq_raw_counts2.txt`: 182,294,665 bytes; MD5 `02a141ff535d351963fd4e8baa3a8629`.\n- Frozen NicheNet ligand-target prior MD5: `b09606b04b2d4490418d9028c0e58b9f`.\n\n## Directory structure\n- `scripts/` — analysis, validation and figure-generation scripts.\n- `results/` — frozen tabular outputs used by the manuscript and figures.\n- `reproducibility/` — console logs, session information and audits.\n- `source_data/` — Nature Communications Source Data workbook, when present.\n- `provenance/` — selected retrieval/model metadata.\n- `MANIFEST_SHA256.csv` — SHA-256 checksum and size for every archived file.\n\n## Key frozen validation results\nIndependent broad-universe NicheNet ranks were concordant between GSE81608 and GSE221156 (Spearman rho = 0.548, P = 1.71e-4). Expression-matched empirical-null sensitivity testing used 10,000 matched random receiver programs per cohort; no ligand achieved empirical FDR < 0.05 in both cohorts. These results are interpreted as sensitivity analyses and regulatory compatibility, not direct signaling evidence.\n\n## Reproducibility notes\nScripts are numbered by workflow stage. Preserve the frozen thresholds, candidate universe, receiver definitions and NicheNet prior described in the manuscript. Do not substitute the GSE81608 1,600-cell representation for the primary 1,492-cell analysis.\n\n## Software environment\nCore analyses were performed in R 4.5.2 on arm64 macOS. Key packages include edgeR 4.8.2, limma 3.66.0, nichenetr 2.2.0, ggplot2 4.0.1 and patchwork 1.3.2. Python 3.14.3 was used for selected download/audit tasks. Exact session records are provided under `reproducibility/` where generated.\n\n## License\nCode in `scripts/` is released under the MIT License. Public-source data and derived outputs remain subject to the terms of their originating repositories/publications.\n'''
(GITHUB / "README.md").write_text(readme, encoding="utf-8")

citation = '''cff-version: 1.2.0\nmessage: "If you use this analysis code, please cite the archived software release and associated manuscript."\ntitle: "Integrated human transcriptomics prioritizes candidate sensory neuron–islet signaling axes in type 2 diabetes: analysis code and reproducibility archive"\ntype: software\nauthors:\n  - family-names: Yadav\n    given-names: Shambhu\n    orcid: "https://orcid.org/0000-0002-0622-1854"\n    affiliation: "Department of Biotechnology, Manipal Institute of Technology, Manipal Academy of Higher Education, Manipal, Karnataka, India"\nversion: 1.0.0\ndate-released: 2026-10-24\nlicense: MIT\nkeywords:\n  - redox biology\n  - sensory neurons\n  - pancreatic beta cells\n  - type 2 diabetes\n  - NicheNet\n  - transcriptomics\n'''
(GITHUB / "CITATION.cff").write_text(citation, encoding="utf-8")

license_text = '''MIT License\n\nCopyright (c) 2026 Shambhu Yadav\n\nPermission is hereby granted, free of charge, to any person obtaining a copy\nof this software and associated documentation files (the "Software"), to deal\nin the Software without restriction, including without limitation the rights\nto use, copy, modify, merge, publish, distribute, sublicense, and/or sell\ncopies of the Software, and to permit persons to whom the Software is\nfurnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\nOUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\nSOFTWARE.\n'''
(GITHUB / "LICENSE").write_text(license_text, encoding="utf-8")

# Ignore local/raw/large files if this folder becomes the repository root.
gitignore = '''# Raw public datasets and large local objects\ndata/raw/\n*.h5ad\n*.h5\n*.rds\n*.RData\n*.Rhistory\n.DS_Store\n.venv/\n__pycache__/\nrelease_v1.0.0/\n'''
(GITHUB / ".gitignore").write_text(gitignore, encoding="utf-8")

# ---------- GitHub manifest ----------
def write_manifest(root, outfile):
    rows=[]
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.name != outfile.name:
            rows.append([str(p.relative_to(root)), p.stat().st_size, sha256(p)])
    with open(outfile, "w", newline="", encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["relative_path","bytes","sha256"]); w.writerows(rows)
    return len(rows)

github_manifest_n = write_manifest(GITHUB, GITHUB / "MANIFEST_SHA256.csv")

# ---------- Zenodo archive ----------
# Start from the GitHub release snapshot.
shutil.copytree(GITHUB, ZENODO / "repository_snapshot", dirs_exist_ok=True)

# Add small processed inputs required for reproduction when available.
processed_count = copy_tree_filtered(
    ROOT / "data/processed", ZENODO / "processed_inputs",
    allowed_ext={".csv", ".tsv", ".txt", ".json", ".rds"}, max_bytes=100*1024*1024
)

# Add final manuscript if present locally (optional, not required for code reproducibility).
for name in [
    "Yadav_Paper3_Nature_Communications_Steps1-3_Integrated.docx",
    "Yadav_Paper3_Nature_Communications_Clean_Manuscript.docx",
]:
    p = ROOT / name
    if p.exists():
        copy_file(p, ZENODO / "manuscript_snapshot" / name)

zenodo_manifest_n = write_manifest(ZENODO, ZENODO / "MANIFEST_SHA256.csv")

# ---------- package-level audit ----------
def count_files(root): return sum(1 for p in root.rglob("*") if p.is_file())
def total_bytes(root): return sum(p.stat().st_size for p in root.rglob("*") if p.is_file())

# Check for files too large for ordinary GitHub repository contents.
large_github=[]
for p in GITHUB.rglob("*"):
    if p.is_file() and p.stat().st_size >= 100*1024*1024:
        large_github.append((str(p.relative_to(GITHUB)), p.stat().st_size))

# Scan text-like files for obvious token/private-key patterns. This is a guardrail, not a guarantee.
secret_patterns = ["ghp_", "github_pat_", "BEGIN PRIVATE KEY", "AWS_SECRET_ACCESS_KEY", "ZENODO_TOKEN"]
secret_hits=[]
for p in GITHUB.rglob("*"):
    if p.is_file() and p.suffix.lower() in {".r",".py",".sh",".txt",".md",".json",".csv",".yml",".yaml",".cff"} and p.stat().st_size < 5*1024*1024:
        try: txt=p.read_text(encoding="utf-8", errors="ignore")
        except Exception: continue
        for pat in secret_patterns:
            if pat in txt:
                secret_hits.append((str(p.relative_to(GITHUB)), pat))

# Build ZIPs.
github_zip = OUT / "Yadav_Paper3_GitHub_Release_v1.0.0.zip"
zenodo_zip = OUT / "Yadav_Paper3_Zenodo_Reproducibility_v1.0.0.zip"
make_zip(GITHUB, github_zip)
make_zip(ZENODO, zenodo_zip)

summary = {
    "release_version":"v1.0.0",
    "scripts_copied":script_count,
    "reproducibility_files_copied":repro_count,
    "result_tables_copied":result_count,
    "processed_inputs_copied_to_zenodo":processed_count,
    "source_data_workbook_included":bool(workbook),
    "github_files":count_files(GITHUB),
    "github_bytes":total_bytes(GITHUB),
    "zenodo_files":count_files(ZENODO),
    "zenodo_bytes":total_bytes(ZENODO),
    "github_files_ge_100MB":large_github,
    "potential_secret_hits":secret_hits,
    "github_zip":str(github_zip),
    "zenodo_zip":str(zenodo_zip),
    "github_zip_sha256":sha256(github_zip),
    "zenodo_zip_sha256":sha256(zenodo_zip),
}
(OUT / "release_audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

print("\n" + "="*72)
print("RELEASE AUDIT")
print("="*72)
print("Scripts copied:", script_count)
print("Reproducibility records copied:", repro_count)
print("Frozen result tables copied:", result_count)
print("Processed inputs added to Zenodo package:", processed_count)
print("Source Data workbook included:", bool(workbook))
print("GitHub package files:", summary["github_files"])
print("GitHub package size: %.2f MB" % (summary["github_bytes"]/1024/1024))
print("Zenodo package files:", summary["zenodo_files"])
print("Zenodo package size: %.2f MB" % (summary["zenodo_bytes"]/1024/1024))
print("GitHub files >=100 MB:", len(large_github))
print("Potential secret-pattern hits:", len(secret_hits))
if secret_hits:
    for x in secret_hits[:20]: print("  REVIEW:", x)
print("\nGitHub ZIP:", github_zip)
print("SHA256:", summary["github_zip_sha256"])
print("\nZenodo ZIP:", zenodo_zip)
print("SHA256:", summary["zenodo_zip_sha256"])
print("\nAudit JSON:", OUT / "release_audit.json")
print("="*72)

if large_github or secret_hits or not workbook:
    print("RELEASE STATUS: REVIEW REQUIRED")
    if not workbook: print("- Source Data workbook missing from release.")
    if large_github: print("- One or more GitHub files are >=100 MB.")
    if secret_hits: print("- Potential secret/token strings require review before upload.")
else:
    print("RELEASE STATUS: PASS")
print("="*72)
