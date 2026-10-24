import json
import os
import sys
import urllib.request

COLLECTION = "58e85c2f-d52e-4c19-8393-b854b84d516e"

API = (
    "https://api.cellxgene.cziscience.com/curation/v1/"
    f"collections/{COLLECTION}"
)

OUTDIR = "data/raw/GSE221156"
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("GSE221156 / CELLxGENE BETA-CELL DATASET DOWNLOAD")
print("=" * 70)

print("\nQuerying CELLxGENE collection...")
with urllib.request.urlopen(API) as r:
    collection = json.load(r)

print("Collection:")
print(collection.get("name", collection.get("title", COLLECTION)))

datasets = collection.get("datasets", [])

print("\nDatasets in collection:", len(datasets))

for i, d in enumerate(datasets, 1):
    title = d.get("title", d.get("name", ""))
    print(
        f"{i}. {title}\n"
        f"   dataset_id: {d.get('dataset_id')}\n"
        f"   cell_count: {d.get('cell_count')}"
    )

# Identify beta-cell-specific dataset
beta = []

for d in datasets:
    title = str(d.get("title", d.get("name", ""))).lower()

    if "beta" in title or "β" in title:
        beta.append(d)

if len(beta) == 0:
    sys.exit(
        "\nERROR: Could not automatically identify beta-cell dataset."
    )

# Prefer dataset closest to expected ~99,029 cells
def cell_distance(d):
    try:
        return abs(int(d.get("cell_count", 0)) - 99029)
    except Exception:
        return 10**12

beta = sorted(beta, key=cell_distance)
d = beta[0]

title = d.get("title", d.get("name", "beta dataset"))
dataset_id = d.get("dataset_id")
cell_count = d.get("cell_count")

print("\nSelected beta-cell dataset:")
print("Title:", title)
print("Dataset ID:", dataset_id)
print("Cell count:", cell_count)

assets = d.get("assets", [])

h5ad = None

for a in assets:
    if str(a.get("filetype", "")).upper() == "H5AD":
        h5ad = a
        break

if h5ad is None:
    sys.exit(
        "\nERROR: No H5AD asset found in collection metadata."
    )

download_url = h5ad.get("url")

if not download_url:
    sys.exit("\nERROR: H5AD asset has no download URL.")

outfile = os.path.join(
    OUTDIR,
    "GSE221156_Beta_processed_CELLxGENE.h5ad"
)

print("\nH5AD asset found.")
print("Output:", outfile)

if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
    print("\nFile already exists; download skipped.")
else:
    print("\nDownloading...")
    urllib.request.urlretrieve(download_url, outfile)
    print("Download complete.")

size = os.path.getsize(outfile)

print("\nDownloaded bytes:", size)
print("Downloaded MB:", round(size / 1024**2, 2))

record = {
    "collection_id": COLLECTION,
    "dataset_id": dataset_id,
    "title": title,
    "cell_count_reported": cell_count,
    "asset_url": download_url,
    "local_file": outfile,
    "size_bytes": size
}

record_file = os.path.join(
    OUTDIR,
    "GSE221156_CELLxGENE_download_record.json"
)

with open(record_file, "w") as f:
    json.dump(record, f, indent=2)

print("\nProvenance record:", record_file)

print("\n" + "=" * 70)
print("SCRIPT 13a COMPLETE")
print("=" * 70)
