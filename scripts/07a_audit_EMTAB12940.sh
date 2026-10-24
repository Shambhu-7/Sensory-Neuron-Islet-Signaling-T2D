#!/bin/bash

set -u

ACC="E-MTAB-12940"
OUT="data/raw/E-MTAB-12940"
REP="reproducibility/07a_EMTAB12940_audit.txt"

mkdir -p "$OUT"
mkdir -p reproducibility

exec > >(tee "$REP") 2>&1

echo "============================================"
echo "SCRIPT 07a: E-MTAB-12940 DATA AUDIT"
echo "============================================"
echo
echo "Accession: $ACC"
echo "Date: $(date)"
echo

INFO_URL="https://www.ebi.ac.uk/biostudies/api/v1/studies/${ACC}/info"
STUDY_URL="https://www.ebi.ac.uk/biostudies/api/v1/studies/${ACC}"

echo "============================================"
echo "1. BIOSTUDIES STUDY INFO"
echo "============================================"

curl -L --fail --retry 5 \
  "$INFO_URL" \
  -o "$OUT/${ACC}_info.json"

if [ $? -ne 0 ]; then
    echo "ERROR: BioStudies info request failed."
    exit 1
fi

python3 - <<'PY'
import json

f = "data/raw/E-MTAB-12940/E-MTAB-12940_info.json"

with open(f) as handle:
    x = json.load(handle)

print("Accession:", x.get("accno", x.get("accession", "NA")))
print("Relative path:", x.get("relPath", "NA"))
print("FTP link:", x.get("ftpLink", "NA"))
print("Release status:", x.get("released", "NA"))
print("File count:", x.get("files", x.get("fileCount", "NA")))
PY

echo
echo "============================================"
echo "2. DOWNLOAD STUDY METADATA JSON"
echo "============================================"

curl -L --fail --retry 5 \
  "$STUDY_URL" \
  -o "$OUT/${ACC}_study.json"

if [ $? -ne 0 ]; then
    echo "WARNING: Detailed study JSON request failed."
fi

echo
echo "============================================"
echo "3. EXTRACT FILE-LIKE RECORDS FROM STUDY JSON"
echo "============================================"

python3 - <<'PY'
import json
import re

f = "data/raw/E-MTAB-12940/E-MTAB-12940_study.json"

try:
    with open(f) as handle:
        obj = json.load(handle)
except Exception as e:
    print("Could not parse study JSON:", e)
    raise SystemExit

hits = set()

extensions = (
    ".rds", ".rda", ".rdata",
    ".csv", ".tsv", ".txt",
    ".mtx", ".h5", ".h5ad",
    ".loom", ".gz", ".zip",
    ".fastq", ".fq"
)

def walk(x):
    if isinstance(x, dict):
        for k, v in x.items():
            walk(k)
            walk(v)
    elif isinstance(x, list):
        for v in x:
            walk(v)
    elif isinstance(x, str):
        s = x.strip()

        low = s.lower()

        if (
            any(ext in low for ext in extensions)
            or "processed" in low
            or "raw" in low
            or "sdrf" in low
            or "idf" in low
        ):
            hits.add(s)

walk(obj)

for x in sorted(hits):
    print(x)
PY

echo
echo "============================================"
echo "4. FTP FILE INVENTORY"
echo "============================================"

FTP_LINK=$(python3 - <<'PY'
import json

f = "data/raw/E-MTAB-12940/E-MTAB-12940_info.json"

with open(f) as handle:
    x = json.load(handle)

print(x.get("ftpLink", ""))
PY
)

echo "FTP base:"
echo "$FTP_LINK"
echo

if [ -n "$FTP_LINK" ]; then

    FTP_FILES="${FTP_LINK%/}/Files/"

    echo "Listing:"
    echo "$FTP_FILES"
    echo

    curl --silent \
         --show-error \
         --list-only \
         "$FTP_FILES" \
         > "$OUT/${ACC}_ftp_file_list.txt"

    echo "Number of FTP entries:"
    wc -l "$OUT/${ACC}_ftp_file_list.txt"

    echo
    echo "First 50 entries:"
    head -50 "$OUT/${ACC}_ftp_file_list.txt"

else
    echo "WARNING: No FTP link found."
fi

echo
echo "============================================"
echo "5. PROCESSED-DATA CANDIDATES"
echo "============================================"

if [ -f "$OUT/${ACC}_ftp_file_list.txt" ]; then

    grep -Ei \
      'rds|rda|rdata|csv|tsv|txt|mtx|h5|h5ad|loom|processed|count|matrix|metadata|annotation|sdrf|idf' \
      "$OUT/${ACC}_ftp_file_list.txt" \
      > "$OUT/${ACC}_processed_candidates.txt" || true

    cat "$OUT/${ACC}_processed_candidates.txt"

else
    echo "No FTP file list available."
fi

echo
echo "============================================"
echo "6. RAW-SEQUENCING INDICATORS"
echo "============================================"

if [ -f "$OUT/${ACC}_ftp_file_list.txt" ]; then

    grep -Ei \
      'fastq|fq.gz|bam|cram' \
      "$OUT/${ACC}_ftp_file_list.txt" \
      | head -30 || true

fi

echo
echo "============================================"
echo "7. DOWNLOAD AUTHORS' ANALYSIS WORKFLOW ONLY"
echo "============================================"

curl -L --fail --retry 5 \
  "https://raw.githubusercontent.com/ManuelMastel/Trace-n-Seq/main/Trace-n-Seq%20SMART-seq2%20R%20Analysis.rmd" \
  -o "$OUT/Trace-n-Seq_SMARTseq2_analysis.rmd"

echo "Workflow file:"
ls -lh "$OUT/Trace-n-Seq_SMARTseq2_analysis.rmd"

echo
echo "Relevant object/file references in workflow:"
grep -Ein \
  'readRDS|saveRDS|raw_SMART|counts_matrix|metadata|pancreas|DRG|JNG' \
  "$OUT/Trace-n-Seq_SMARTseq2_analysis.rmd" \
  | head -80

echo
echo "============================================"
echo "SCRIPT 07a COMPLETED"
echo "============================================"
echo
echo "NO LARGE EXPRESSION OR FASTQ DATA WERE"
echo "INTENTIONALLY DOWNLOADED BY THIS AUDIT."
echo "============================================"

