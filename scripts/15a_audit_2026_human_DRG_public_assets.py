#!/usr/bin/env python3
"""
15a_audit_2026_human_DRG_public_assets.py

Purpose
-------
Audit public availability of the 2026 human DRG diabetes/DPN multi-omic cohort:
"Progressive neurodegeneration in human dorsal root ganglion from diabetes to painful neuropathy"
bioRxiv DOI: 10.64898/2026.01.16.700028
PMC: PMC12871725

This script DOES NOT download raw sequencing data or large processed matrices.
It only queries public metadata/APIs and saves small JSON/XML/text responses.

Known SPARC DOIs from the paper:
  10.26275/vz0c-uet7  -> previously published Visium dataset (GSE295206 lineage)
  10.26275/0jnf-7wq9  -> candidate second public dataset; availability to be audited

No third-party Python packages are required.
"""

from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError, URLError
import json
import re
import ssl
import time
import hashlib
import sys

ROOT = Path(".")
OUT = ROOT / "results" / "validation" / "DRG2026_public_asset_audit"
RAW = ROOT / "data" / "raw" / "DRG2026_public_asset_audit"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

UA = "Mozilla/5.0 Paper3-Sensory-Islet reproducibility audit; contact=local-research-script"
CTX = ssl.create_default_context()

ARTICLE_DOI = "10.64898/2026.01.16.700028"
PMC_ID = "PMC12871725"
SPARC_DOIS = ["10.26275/vz0c-uet7", "10.26275/0jnf-7wq9"]

def fetch(url, timeout=45, accept="application/json"):
    req = Request(url, headers={"User-Agent": UA, "Accept": accept})
    try:
        with urlopen(req, timeout=timeout, context=CTX) as r:
            body = r.read()
            return {
                "ok": True,
                "status": getattr(r, "status", 200),
                "final_url": r.geturl(),
                "content_type": r.headers.get("Content-Type", ""),
                "body": body,
            }
    except HTTPError as e:
        try:
            body = e.read()
        except Exception:
            body = b""
        return {
            "ok": False, "status": e.code, "final_url": url,
            "content_type": e.headers.get("Content-Type", "") if e.headers else "",
            "body": body, "error": str(e)
        }
    except Exception as e:
        return {
            "ok": False, "status": None, "final_url": url,
            "content_type": "", "body": b"", "error": repr(e)
        }

def save_bytes(path, body):
    path.write_bytes(body)

def safe_json(body):
    try:
        return json.loads(body.decode("utf-8", errors="replace"))
    except Exception:
        return None

def md5_bytes(body):
    return hashlib.md5(body).hexdigest()

def recursive_find(obj, keys):
    """Find dictionaries containing any requested key; useful for unknown API response wrappers."""
    hits = []
    if isinstance(obj, dict):
        if any(k in obj for k in keys):
            hits.append(obj)
        for v in obj.values():
            hits.extend(recursive_find(v, keys))
    elif isinstance(obj, list):
        for v in obj:
            hits.extend(recursive_find(v, keys))
    return hits

def compact_dataset_records(obj):
    hits = recursive_find(obj, {"id", "name", "doi", "version", "versionId"})
    records = []
    seen = set()
    for h in hits:
        did = h.get("id", h.get("datasetId"))
        name = h.get("name", h.get("datasetName", h.get("title")))
        doi = h.get("doi")
        version = h.get("version", h.get("versionId"))
        if did is None and name is None and doi is None:
            continue
        key = (str(did), str(name), str(doi), str(version))
        if key in seen:
            continue
        seen.add(key)
        records.append({
            "id": did, "name": name, "doi": doi, "version": version,
            "description": h.get("description"),
            "size": h.get("size"),
            "status": h.get("status"),
        })
    return records

print("=" * 78)
print("SCRIPT 15a — 2026 HUMAN DRG PUBLIC-ASSET AUDIT")
print("=" * 78)
print("Article DOI:", ARTICLE_DOI)
print("PMC:", PMC_ID)
print()
print("Prespecified purpose:")
print("  Orthogonal disease-context validation of the Paper 3 DRG findings.")
print("  NOT an independent replication claim for GSE295206.")
print("  NO raw sequencing or large matrix downloads in this script.")
print()

# -------------------------------------------------------------------------
# 1. DataCite metadata for the two SPARC DOIs
# -------------------------------------------------------------------------
print("1) DATACITE DOI AUDIT")
doi_rows = []
for doi in SPARC_DOIS:
    url = "https://api.datacite.org/dois/" + doi
    res = fetch(url)
    outpath = RAW / ("datacite_" + doi.replace("/", "_").replace(".", "_") + ".json")
    save_bytes(outpath, res["body"])
    j = safe_json(res["body"])
    title = None
    publisher = None
    state = None
    landing = None
    if isinstance(j, dict):
        attrs = (j.get("data") or {}).get("attributes") or {}
        titles = attrs.get("titles") or []
        if titles and isinstance(titles[0], dict):
            title = titles[0].get("title")
        publisher = attrs.get("publisher")
        state = attrs.get("state")
        landing = attrs.get("url")
    doi_rows.append({
        "doi": doi, "http_status": res["status"], "title": title,
        "publisher": publisher, "state": state, "landing_url": landing
    })
    print(f"  {doi}")
    print(f"    HTTP: {res['status']}")
    print(f"    title: {title}")
    print(f"    publisher: {publisher}")
    print(f"    state: {state}")
    print(f"    landing: {landing}")

(OUT / "SPARC_DOI_audit.json").write_text(json.dumps(doi_rows, indent=2), encoding="utf-8")
print()

# -------------------------------------------------------------------------
# 2. DOI resolver audit
# -------------------------------------------------------------------------
print("2) DOI RESOLUTION AUDIT")
resolver_rows = []
for doi in SPARC_DOIS:
    url = "https://doi.org/" + doi
    res = fetch(url, accept="text/html,application/xhtml+xml")
    resolver_rows.append({
        "doi": doi, "http_status": res["status"], "final_url": res["final_url"],
        "content_type": res["content_type"]
    })
    print(f"  {doi}")
    print(f"    HTTP: {res['status']}")
    print(f"    final URL: {res['final_url']}")
(OUT / "SPARC_DOI_resolution.json").write_text(json.dumps(resolver_rows, indent=2), encoding="utf-8")
print()

# -------------------------------------------------------------------------
# 3. Pennsieve Discover searches
# -------------------------------------------------------------------------
print("3) PENNSIEVE DISCOVER SEARCH")
search_terms = [
    "Progressive neurodegeneration human dorsal root ganglion",
    "dorsal root ganglion diabetes DPN",
    "single nucleus dorsal root ganglion diabetes",
    "Sankaranarayanan dorsal root ganglion",
    "Price dorsal root ganglion diabetes",
]
all_records = []
for i, term in enumerate(search_terms, 1):
    params = urlencode({
        "limit": 50,
        "offset": 0,
        "query": term,
        "orderBy": "relevance",
        "orderDirection": "desc"
    })
    url = "https://api.pennsieve.io/discover/search/datasets?" + params
    res = fetch(url)
    save_bytes(RAW / f"pennsieve_search_{i}.json", res["body"])
    j = safe_json(res["body"])
    records = compact_dataset_records(j) if j is not None else []
    print(f"  Query {i}: {term}")
    print(f"    HTTP: {res['status']}; candidate records parsed: {len(records)}")
    for r in records[:10]:
        print(f"      id={r['id']} | version={r['version']} | doi={r['doi']} | {r['name']}")
        all_records.append(r)
    time.sleep(0.2)

# Deduplicate records
dedup = []
seen = set()
for r in all_records:
    key = (str(r.get("id")), str(r.get("name")), str(r.get("doi")))
    if key not in seen:
        seen.add(key)
        dedup.append(r)

(OUT / "Pennsieve_search_candidates.json").write_text(
    json.dumps(dedup, indent=2, default=str), encoding="utf-8"
)
print(f"\n  Unique parsed candidate records: {len(dedup)}")
print()

# -------------------------------------------------------------------------
# 4. Explicitly inspect known dataset 490 (the prior Visium dataset)
# -------------------------------------------------------------------------
print("4) KNOWN GSE295206-LINEAGE SPARC DATASET CONTROL")
known_url = "https://api.pennsieve.io/discover/datasets/490"
res490 = fetch(known_url)
save_bytes(RAW / "pennsieve_dataset_490.json", res490["body"])
j490 = safe_json(res490["body"])
print("  Dataset 490 HTTP:", res490["status"])
if isinstance(j490, dict):
    print("  Dataset 490 response saved for provenance.")
print("  Expected identity: spatial transcriptomics / prior GSE295206 lineage.")
print()

# -------------------------------------------------------------------------
# 5. Try to identify candidate dataset ID for second DOI from search/metadata
# -------------------------------------------------------------------------
print("5) SECOND SPARC DOI IDENTIFICATION")
target = "0jnf-7wq9"
candidate_ids = set()

for r in dedup:
    blob = json.dumps(r, default=str).lower()
    if target in blob or "progressive neurodegeneration" in blob or "single nucleus" in blob:
        if r.get("id") is not None:
            try:
                candidate_ids.add(int(r["id"]))
            except Exception:
                pass

# Also inspect DataCite landing URL for a /datasets/<id> pattern
for row in doi_rows:
    if row["doi"].lower().endswith(target):
        landing = row.get("landing_url") or ""
        m = re.search(r"/datasets/(\d+)", landing)
        if m:
            candidate_ids.add(int(m.group(1)))

if candidate_ids:
    print("  Candidate Pennsieve dataset IDs:", sorted(candidate_ids))
else:
    print("  No Pennsieve dataset ID could yet be identified for 10.26275/0jnf-7wq9.")

# -------------------------------------------------------------------------
# 6. Inspect candidate dataset metadata + root files only
# -------------------------------------------------------------------------
candidate_summaries = []
for did in sorted(candidate_ids):
    meta_url = f"https://api.pennsieve.io/discover/datasets/{did}"
    res = fetch(meta_url)
    save_bytes(RAW / f"pennsieve_dataset_{did}.json", res["body"])
    j = safe_json(res["body"])

    # Try to infer current version
    version = None
    if isinstance(j, dict):
        candidates = recursive_find(j, {"version", "versionId"})
        for h in candidates:
            v = h.get("version", h.get("versionId"))
            try:
                version = int(v)
                break
            except Exception:
                pass
    if version is None:
        version = 1

    browse_url = (
        f"https://api.pennsieve.io/discover/datasets/{did}/versions/{version}/files/browse"
        "?limit=200&offset=0"
    )
    br = fetch(browse_url)
    save_bytes(RAW / f"pennsieve_dataset_{did}_v{version}_root_files.json", br["body"])
    bj = safe_json(br["body"])

    # Extract likely filenames/path strings without assuming exact schema.
    strings = []
    def collect_strings(x):
        if isinstance(x, str):
            strings.append(x)
        elif isinstance(x, dict):
            for vv in x.values():
                collect_strings(vv)
        elif isinstance(x, list):
            for vv in x:
                collect_strings(vv)
    if bj is not None:
        collect_strings(bj)
    likely = sorted(set(
        s for s in strings
        if re.search(r"\.(rds|h5ad|h5|mtx|csv|tsv|txt|xlsx|r|rmd|qs|parquet)(\.gz)?$", s, re.I)
    ))
    candidate_summaries.append({
        "dataset_id": did,
        "metadata_http": res["status"],
        "version_used_for_browse": version,
        "browse_http": br["status"],
        "likely_processed_or_code_files": likely[:500],
    })
    print(f"\n  Dataset {did}: metadata HTTP {res['status']}; root browse v{version} HTTP {br['status']}")
    if likely:
        print("  Likely processed/code files found at root:")
        for s in likely[:40]:
            print("    ", s)
    else:
        print("  No obvious processed matrix/table/code filenames parsed at root.")

(OUT / "Candidate_dataset_file_audit.json").write_text(
    json.dumps(candidate_summaries, indent=2), encoding="utf-8"
)
print()

# -------------------------------------------------------------------------
# 7. PMC Open Access package metadata (metadata only; no tarball download)
# -------------------------------------------------------------------------
print("6) PMC OPEN-ACCESS PACKAGE AUDIT")
oa_url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={PMC_ID}"
oa = fetch(oa_url, accept="application/xml,text/xml")
save_bytes(RAW / f"{PMC_ID}_oa_manifest.xml", oa["body"])
oa_text = oa["body"].decode("utf-8", errors="replace")
print("  OA API HTTP:", oa["status"])
links = re.findall(r'href="([^"]+)"', oa_text)
if links:
    print("  Package links found:")
    for x in links:
        print("   ", x)
else:
    print("  No OA package links parsed.")
print("  Large OA package is NOT downloaded by this script.")
print()

# -------------------------------------------------------------------------
# 8. Sensoryomics reachability
# -------------------------------------------------------------------------
print("7) SENSORYOMICS SUMMARY-DATA SITE AUDIT")
for url in ["https://sensoryomics.com", "https://www.sensoryomics.com"]:
    res = fetch(url, accept="text/html")
    print(f"  {url} -> HTTP {res['status']} -> {res['final_url']}")
    if res["ok"]:
        save_bytes(RAW / "sensoryomics_homepage.html", res["body"])
        break
print()

# -------------------------------------------------------------------------
# 9. Decision summary
# -------------------------------------------------------------------------
second = next((x for x in doi_rows if x["doi"].lower().endswith(target)), None)
second_registered = bool(second and second.get("http_status") == 200 and second.get("title"))

print("=" * 78)
print("AUDIT DECISION")
print("=" * 78)
print("Known Visium DOI 10.26275/vz0c-uet7:")
print("  Treat as GSE295206-lineage data; NOT independent validation.")
print()
print("Second DOI 10.26275/0jnf-7wq9:")
if second_registered:
    print("  DOI metadata is registered/public.")
else:
    print("  DOI metadata was not clearly public/registered through DataCite at audit time.")

if candidate_ids:
    print("  A candidate Pennsieve dataset ID was identified:", ", ".join(map(str, sorted(candidate_ids))))
    print("  Review the printed/root-file audit before downloading anything.")
else:
    print("  No live Pennsieve dataset ID was identified from public search/API metadata.")
    print("  In that case, do NOT infer that the full snRNA-seq matrix is publicly downloadable yet.")
print()
print("Next analysis, ONLY if suitable processed snRNA-seq or summary tables are public:")
print("  1. Prespecified 13 sensory-source ligand expression in human neurons.")
print("  2. Main candidates: SEMA7A, EDIL3, FLRT3, PGF; SPP1 reference.")
print("  3. Prespecified redox anchors/pathway-direction concordance.")
print("  4. Control -> Diabetes -> DPN disease-context behavior.")
print("  5. Donor-aware statistics only when donor-resolved data permit them.")
print()
print("Guardrail:")
print("  Call this 'orthogonal disease-context validation in an expanded human DRG")
print("  single-nucleus cohort', not 'independent replication of GSE295206'.")
print()
print("Saved audit directory:", OUT)
print("Saved raw API responses:", RAW)
print("=" * 78)
print("SCRIPT 15a COMPLETE")
print("=" * 78)
