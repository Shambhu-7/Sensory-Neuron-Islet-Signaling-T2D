#!/usr/bin/env python3
"""
15b_recover_2026_DRG_supplementary_tables.py

Audit PMC/NCBI/Europe PMC article representations for supplementary attachments
from PMC12871725 and download ONLY small tabular/archive supplementary files.

No raw sequencing data are downloaded.
"""

from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urljoin, urlparse
from urllib.error import HTTPError
import re, json, hashlib, ssl, os

PMC = "PMC12871725"
OUT = Path("results/validation/DRG2026_supplementary_audit")
RAW = Path("data/raw/DRG2026_supplementary")
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

UA = "Mozilla/5.0 Paper3-Sensory-Islet supplementary-data audit"
CTX = ssl.create_default_context()
MAX_DOWNLOAD = 75 * 1024 * 1024  # 75 MB
TABULAR_EXT = (".xlsx", ".xls", ".csv", ".tsv", ".txt", ".zip", ".docx")

sources = [
    ("pmc_html", f"https://pmc.ncbi.nlm.nih.gov/articles/{PMC}/"),
    ("ncbi_efetch_xml",
     f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={PMC.replace('PMC','')}&retmode=xml"),
    ("europepmc_xml",
     f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMC}/fullTextXML"),
]

def fetch(url, method="GET", headers=None, timeout=60):
    h = {"User-Agent": UA}
    if headers:
        h.update(headers)
    req = Request(url, headers=h, method=method)
    try:
        with urlopen(req, timeout=timeout, context=CTX) as r:
            return {
                "ok": True,
                "status": getattr(r, "status", 200),
                "url": r.geturl(),
                "headers": dict(r.headers.items()),
                "body": r.read() if method != "HEAD" else b"",
            }
    except HTTPError as e:
        try:
            body = e.read()
        except Exception:
            body = b""
        return {
            "ok": False, "status": e.code, "url": url,
            "headers": dict(e.headers.items()) if e.headers else {},
            "body": body, "error": str(e)
        }
    except Exception as e:
        return {"ok": False, "status": None, "url": url,
                "headers": {}, "body": b"", "error": repr(e)}

def clean_href(x):
    return (x.replace("&amp;", "&")
             .replace("&#x2F;", "/")
             .replace("&#47;", "/")
             .strip())

print("="*78)
print("SCRIPT 15b — RECOVER 2026 HUMAN DRG SUPPLEMENTARY TABLES")
print("="*78)
print("Article:", PMC)
print("Policy: download only small supplementary tables/archives; no raw sequencing.\n")

all_links = []
source_report = []

for label, url in sources:
    print(f"FETCH {label}: {url}")
    r = fetch(url)
    suffix = ".xml" if "xml" in label else ".html"
    path = RAW / f"{label}{suffix}"
    path.write_bytes(r["body"])
    print(f"  HTTP: {r['status']}; bytes: {len(r['body']):,}")
    source_report.append({
        "label": label, "url": url, "status": r["status"],
        "bytes": len(r["body"]), "saved": str(path)
    })
    if not r["body"]:
        continue

    text = r["body"].decode("utf-8", errors="replace")

    # HTML/XML hrefs, including xlink:href.
    hrefs = re.findall(r'(?:href|xlink:href)\s*=\s*["\']([^"\']+)["\']', text, flags=re.I)

    # Also catch bare supplementary-looking file references.
    bare = re.findall(
        r'(?i)([A-Za-z0-9_./:%?=&+\-]+(?:supp(?:lement)?[A-Za-z0-9_./:%?=&+\-]*|media-[0-9]+|nihpp-[A-Za-z0-9_.\-]+)'
        r'\.(?:xlsx|xls|csv|tsv|txt|zip|docx|pdf))',
        text
    )
    hrefs.extend(bare)

    base = r["url"] if r["ok"] else url
    for h in hrefs:
        h = clean_href(h)
        absu = urljoin(base, h)
        low = absu.lower()
        if any(k in low for k in [
            "supp", "media-", "/bin/", ".xlsx", ".xls", ".csv",
            ".tsv", ".zip", ".docx", "nihpp"
        ]):
            all_links.append({"source": label, "url": absu})

# Deduplicate.
uniq = []
seen = set()
for x in all_links:
    u = x["url"]
    if u not in seen:
        seen.add(u)
        uniq.append(x)

print(f"\nCandidate supplementary links discovered: {len(uniq)}")

audit = []
for i, x in enumerate(uniq, 1):
    u = x["url"]
    print(f"\n[{i}] {u}")

    # HEAD first.
    h = fetch(u, method="HEAD", timeout=45)
    status = h["status"]
    headers = h["headers"]
    final_url = h["url"]

    # Some servers reject HEAD: use 1-byte range GET only for metadata.
    if not h["ok"]:
        h = fetch(u, method="GET", headers={"Range": "bytes=0-0"}, timeout=45)
        status = h["status"]
        headers = h["headers"]
        final_url = h["url"]

    ctype = headers.get("Content-Type", headers.get("content-type", ""))
    clen = headers.get("Content-Length", headers.get("content-length"))
    cdisp = headers.get("Content-Disposition", headers.get("content-disposition", ""))

    try:
        size = int(clen) if clen else None
    except Exception:
        size = None

    parsed_name = os.path.basename(urlparse(final_url).path) or f"supplement_{i}"
    m = re.search(r'filename\*?=(?:UTF-8\'\')?["\']?([^"\';]+)', cdisp, flags=re.I)
    if m:
        parsed_name = m.group(1).strip()

    lowname = parsed_name.lower()
    is_tabular = lowname.endswith(TABULAR_EXT)
    # If filename has no useful extension, infer from URL.
    if not is_tabular:
        lowurl = final_url.lower().split("?")[0]
        is_tabular = lowurl.endswith(TABULAR_EXT)

    print(f"    HTTP: {status}")
    print(f"    final: {final_url}")
    print(f"    content-type: {ctype}")
    print(f"    content-length: {size}")
    print(f"    filename: {parsed_name}")
    print(f"    tabular/archive candidate: {is_tabular}")

    downloaded = False
    saved = None
    md5 = None

    if status in (200, 206) and is_tabular and (size is None or size <= MAX_DOWNLOAD):
        # Full GET.
        g = fetch(final_url, method="GET", timeout=180)
        body = g["body"]
        # Guard against accidentally receiving a huge object.
        if g["ok"] and len(body) <= MAX_DOWNLOAD:
            # Sanitize filename.
            parsed_name = re.sub(r'[^A-Za-z0-9._+\-]+', '_', parsed_name)
            if not parsed_name or parsed_name.startswith("_"):
                parsed_name = f"supplement_{i}.bin"
            dest = RAW / parsed_name
            # Avoid overwriting different files with same basename.
            if dest.exists() and dest.read_bytes() != body:
                dest = RAW / f"{i:02d}_{parsed_name}"
            dest.write_bytes(body)
            downloaded = True
            saved = str(dest)
            md5 = hashlib.md5(body).hexdigest()
            print(f"    DOWNLOADED: {dest} ({len(body):,} bytes)")
            print(f"    MD5: {md5}")
        elif len(body) > MAX_DOWNLOAD:
            print("    SKIPPED: response exceeded 75 MB safety limit.")
        else:
            print(f"    Full GET failed: HTTP {g['status']}")
    else:
        if is_tabular and size and size > MAX_DOWNLOAD:
            print("    SKIPPED: >75 MB safety limit.")
        else:
            print("    Not downloaded.")

    audit.append({
        "source": x["source"], "original_url": u, "http_status": status,
        "final_url": final_url, "content_type": ctype, "content_length": size,
        "filename": parsed_name, "tabular_archive_candidate": is_tabular,
        "downloaded": downloaded, "saved": saved, "md5": md5
    })

(OUT / "article_source_fetch_report.json").write_text(
    json.dumps(source_report, indent=2), encoding="utf-8"
)
(OUT / "supplementary_link_audit.json").write_text(
    json.dumps(audit, indent=2), encoding="utf-8"
)

# Inventory downloaded files.
downloaded = [x for x in audit if x["downloaded"]]
print("\n" + "="*78)
print("DOWNLOADED SUPPLEMENTARY FILE INVENTORY")
print("="*78)
if not downloaded:
    print("No downloadable small tabular/archive supplementary files were recovered.")
else:
    for x in downloaded:
        print(f"  {x['saved']}")
        print(f"    bytes={Path(x['saved']).stat().st_size:,} md5={x['md5']}")

# Scan filenames and lightweight contents for likely Supplemental Data numbering.
print("\nLIKELY HIGH-VALUE TARGETS")
print("  Supplementary Data 3: Control vs Diabetes pseudobulk DE")
print("  Supplementary Data 5: neuronal Diabetes vs DPN DE")
print("  Supplementary Data 6: neuronal Nageotte-burden DE")
print("  Supplementary Data 8: non-neuronal ligand DE")
print("  Supplementary Data 9: neuronal receptor DE")
print("  Supplementary Data 11/12: Xenium subtype/DE validation")

print("\nDECISION")
if downloaded:
    print("  PASS: one or more small supplementary assets were recovered.")
    print("  Next: inspect table/sheet names and map them to Supplementary Data numbers.")
else:
    print("  No small supplementary tables were recoverable from the current public article endpoints.")
    print("  Combined with Script 15a, the 2026 snRNA cohort should then be treated as")
    print("  literature-level orthogonal support only until the promised processed data are released.")

print("\nSaved audit:", OUT / "supplementary_link_audit.json")
print("="*78)
print("SCRIPT 15b COMPLETE")
print("="*78)
