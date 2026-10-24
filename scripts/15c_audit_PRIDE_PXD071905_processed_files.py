#!/usr/bin/env python3
"""
15c_audit_PRIDE_PXD071905_processed_files.py

Audit PRIDE/ProteomeXchange project PXD071905 for small processed result files
from the 2026 human DRG diabetes/DPN study.

This script does NOT download raw MS files, mzML files, or large archives.
It only queries public metadata/file listings and writes small JSON records.
No third-party Python packages required.
"""

from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError
import json, ssl, re, hashlib

ACC = "PXD071905"
OUT = Path("results/validation/DRG2026_PRIDE_audit")
RAW = Path("data/raw/DRG2026_PRIDE_audit")
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

UA = "Mozilla/5.0 Paper3-Sensory-Islet PRIDE processed-file audit"
CTX = ssl.create_default_context()

def fetch(url, timeout=60):
    req = Request(url, headers={"User-Agent": UA, "Accept": "application/json,text/plain,*/*"})
    try:
        with urlopen(req, timeout=timeout, context=CTX) as r:
            return {
                "ok": True, "status": getattr(r, "status", 200),
                "final_url": r.geturl(), "headers": dict(r.headers.items()),
                "body": r.read()
            }
    except HTTPError as e:
        try:
            body = e.read()
        except Exception:
            body = b""
        return {"ok": False, "status": e.code, "final_url": url,
                "headers": dict(e.headers.items()) if e.headers else {},
                "body": body, "error": str(e)}
    except Exception as e:
        return {"ok": False, "status": None, "final_url": url,
                "headers": {}, "body": b"", "error": repr(e)}

def safe_json(body):
    try:
        return json.loads(body.decode("utf-8", errors="replace"))
    except Exception:
        return None

def collect_filelike(obj, out):
    if isinstance(obj, dict):
        keys = {str(k).lower() for k in obj.keys()}
        if keys & {"filename","file_name","name","downloadlink","download_link","filetype","file_type"}:
            out.append(obj)
        for v in obj.values():
            collect_filelike(v, out)
    elif isinstance(obj, list):
        for v in obj:
            collect_filelike(v, out)

def first(d, names):
    for n in names:
        if n in d and d[n] not in (None, ""):
            return d[n]
    # case-insensitive
    low = {str(k).lower(): k for k in d}
    for n in names:
        if n.lower() in low:
            return d[low[n.lower()]]
    return None

print("="*78)
print("SCRIPT 15c — PRIDE PXD071905 PROCESSED-FILE AUDIT")
print("="*78)
print("Accession:", ACC)
print("Purpose: orthogonal human DRG proteomic disease-context validation.")
print("Safety: metadata/file-list audit only; NO raw MS download.\n")

endpoints = [
    ("project_v3", f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACC}"),
    ("files_v3", f"https://www.ebi.ac.uk/pride/ws/archive/v3/files/byProject?accession={ACC}"),
    ("project_v2", f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{ACC}"),
    ("files_v2", f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={ACC}"),
]

responses = {}
for label, url in endpoints:
    r = fetch(url)
    responses[label] = r
    ext = ".json"
    (RAW / f"{label}{ext}").write_bytes(r["body"])
    print(f"{label}: HTTP {r['status']} | bytes {len(r['body']):,} | {r['final_url']}")

print()

# Parse project metadata from whichever endpoint works.
project = None
for label in ("project_v3","project_v2"):
    j = safe_json(responses[label]["body"])
    if isinstance(j, dict):
        project = j
        break

if project:
    title = first(project, ["title","projectTitle"])
    submission = first(project, ["submissionDate","submission_date"])
    publication = first(project, ["publicationDate","publication_date"])
    print("PROJECT METADATA")
    print("  title:", title)
    print("  submission:", submission)
    print("  publication:", publication)
else:
    print("PROJECT METADATA: no parseable public project JSON returned.")

# Parse files.
raw_records = []
for label in ("files_v3","files_v2"):
    j = safe_json(responses[label]["body"])
    if j is not None:
        collect_filelike(j, raw_records)

# Deduplicate and normalize.
files = []
seen = set()
for d in raw_records:
    name = first(d, ["fileName","filename","file_name","name"])
    if not name:
        continue
    name = str(name)
    size = first(d, ["fileSizeBytes","fileSize","size","file_size"])
    ftype = first(d, ["fileType","filetype","file_type"])
    dl = first(d, ["downloadLink","download_link","downloadUrl","downloadURL","ftpLink","publicFileLocations"])
    key = name
    if key in seen:
        continue
    seen.add(key)
    files.append({"name": name, "size": size, "type": ftype, "download": dl})

print("\nFILES PARSED:", len(files))

def category(name):
    x = name.lower()
    raw_ext = (".raw",".d",".wiff",".wiff.scan",".mzml",".mzxml")
    if x.endswith(raw_ext) or "raw" in x and not any(t in x for t in ["report","result","protein","matrix"]):
        return "RAW_MS"
    if any(k in x for k in [
        "report", "protein", "proteingroup", "protein_group", "matrix",
        "diann", "dia-nn", "pg_matrix", "pr_matrix", "quant",
        "result", "output", "statistics", "supp", "metadata"
    ]) or x.endswith((".tsv",".csv",".txt",".xlsx",".xls",".parquet")):
        return "PROCESSED_CANDIDATE"
    if x.endswith((".zip",".gz",".tar",".tgz")):
        return "ARCHIVE_UNKNOWN"
    return "OTHER"

for f in files:
    f["category"] = category(f["name"])

counts = {}
for f in files:
    counts[f["category"]] = counts.get(f["category"], 0) + 1

print("FILE CATEGORIES")
for k in sorted(counts):
    print(f"  {k}: {counts[k]}")

processed = [f for f in files if f["category"] == "PROCESSED_CANDIDATE"]
print("\nPROCESSED CANDIDATES")
if processed:
    for i, f in enumerate(processed, 1):
        print(f"  {i}. {f['name']}")
        print(f"     size={f['size']} type={f['type']}")
        print(f"     download={f['download']}")
else:
    print("  None identified from public PRIDE file listing.")

(OUT / "PXD071905_file_inventory.json").write_text(
    json.dumps(files, indent=2, default=str), encoding="utf-8"
)
(OUT / "PXD071905_processed_candidates.json").write_text(
    json.dumps(processed, indent=2, default=str), encoding="utf-8"
)

print("\nDECISION")
if processed:
    print("  PASS: processed candidate files are publicly listed.")
    print("  Next step: inspect/download ONLY the smallest useful DIA-NN/protein table(s).")
    print("  Do NOT download RAW/mzML files.")
else:
    print("  No obvious processed protein table is currently exposed through the public PRIDE API.")
    print("  Do not download raw MS files merely to add this validation.")

print("\nGuardrail:")
print("  This cohort can support orthogonal disease-context/proteomic validation,")
print("  not independent replication of the GSE295206 spatial cohort.")

print("\nSaved:")
print(" ", OUT / "PXD071905_file_inventory.json")
print(" ", OUT / "PXD071905_processed_candidates.json")
print("="*78)
print("SCRIPT 15c COMPLETE")
print("="*78)
