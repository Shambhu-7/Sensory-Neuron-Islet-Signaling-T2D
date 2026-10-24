#!/usr/bin/env python3
import os, sys, json, math

FILE = "data/raw/GSE221156/GSE221156_Beta_processed_CELLxGENE.h5ad"
OUT = "results/validation/GSE221156/GSE221156_H5AD_structure_audit.txt"

try:
    import h5py
    import numpy as np
except ImportError as e:
    sys.exit(
        "\nERROR: This audit needs h5py and numpy.\n"
        "Install once with:\n"
        "  python3 -m pip install --user h5py numpy\n"
        "Then rerun this script.\n"
    )

os.makedirs(os.path.dirname(OUT), exist_ok=True)

def dec(x):
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="replace")
    if isinstance(x, np.bytes_):
        return x.astype(str)
    return str(x)

def attr_dict(obj):
    out = {}
    for k, v in obj.attrs.items():
        try:
            if isinstance(v, np.ndarray):
                out[k] = v.tolist()
            elif isinstance(v, (bytes, np.bytes_)):
                out[k] = dec(v)
            else:
                out[k] = v
        except Exception:
            out[k] = str(v)
    return out

def shape_of(obj):
    if isinstance(obj, h5py.Dataset):
        return tuple(obj.shape)
    if isinstance(obj, h5py.Group):
        if "shape" in obj.attrs:
            try:
                return tuple(int(x) for x in obj.attrs["shape"])
            except Exception:
                pass
        # sparse matrix fallback
        if "indptr" in obj and "indices" in obj:
            return "sparse; shape attr unavailable"
    return None

def sample_numeric_dataset(ds, n=200000):
    if not isinstance(ds, h5py.Dataset) or ds.size == 0:
        return None
    m = min(int(ds.size), n)
    arr = ds[:m]
    arr = np.asarray(arr)
    if arr.dtype.kind not in "iuf":
        return None
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    nz = arr[arr != 0]
    z = {
        "sample_n": int(arr.size),
        "sample_min": float(np.min(arr)),
        "sample_max": float(np.max(arr)),
        "sample_nonzero_n": int(nz.size),
    }
    if nz.size:
        frac = np.abs(nz - np.round(nz))
        z["sample_nonzero_integer_like_fraction"] = float(np.mean(frac < 1e-8))
        qs = np.quantile(nz, [0, .25, .5, .75, .9, .99, 1])
        z["sample_nonzero_quantiles"] = [float(x) for x in qs]
    return z

def matrix_info(f, path):
    if path not in f:
        return None
    obj = f[path]
    info = {
        "path": path,
        "type": "dataset" if isinstance(obj, h5py.Dataset) else "group",
        "shape": shape_of(obj),
        "attrs": attr_dict(obj)
    }
    if isinstance(obj, h5py.Dataset):
        info["dtype"] = str(obj.dtype)
        info["numeric_sample"] = sample_numeric_dataset(obj)
    elif isinstance(obj, h5py.Group):
        info["keys"] = list(obj.keys())
        if "data" in obj:
            info["data_dtype"] = str(obj["data"].dtype)
            info["numeric_sample"] = sample_numeric_dataset(obj["data"])
            info["stored_values"] = int(obj["data"].size)
    return info

def read_vector(g, key, max_full=500000):
    if key not in g:
        return None
    obj = g[key]

    # AnnData categorical
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        codes = np.asarray(obj["codes"][:], dtype=int)
        cats = [dec(x) for x in obj["categories"][:]]
        vals = np.array(
            [cats[c] if 0 <= c < len(cats) else "<NA>" for c in codes],
            dtype=object
        )
        return vals

    # Nullable/string-like group
    if isinstance(obj, h5py.Group):
        if "values" in obj:
            vals = obj["values"][:]
            return np.array([dec(x) for x in vals], dtype=object)
        return None

    if isinstance(obj, h5py.Dataset):
        if obj.size > max_full:
            return None
        vals = obj[:]
        if vals.dtype.kind in ("S", "O", "U"):
            return np.array([dec(x) for x in vals], dtype=object)
        return np.asarray(vals)

    return None

def value_counts(vals, top=60):
    if vals is None:
        return None
    vals2 = []
    for x in vals:
        if x is None:
            vals2.append("<NA>")
        elif isinstance(x, float) and np.isnan(x):
            vals2.append("<NA>")
        else:
            vals2.append(str(x))
    u, c = np.unique(np.asarray(vals2, dtype=object), return_counts=True)
    order = np.argsort(-c)
    return [(str(u[i]), int(c[i])) for i in order[:top]]

lines = []
def p(*args):
    s = " ".join(str(x) for x in args)
    print(s)
    lines.append(s)

p("="*78)
p("GSE221156 / CELLxGENE H5AD STRUCTURE + ASSAY AUDIT")
p("="*78)
p("File:", FILE)
p("Bytes:", os.path.getsize(FILE))

with h5py.File(FILE, "r") as f:
    p("\nTOP-LEVEL KEYS")
    p(list(f.keys()))

    p("\nROOT ATTRIBUTES")
    p(json.dumps(attr_dict(f), indent=2, default=str))

    for path in ["X", "raw/X"]:
        info = matrix_info(f, path)
        p(f"\nMATRIX AUDIT: {path}")
        if info is None:
            p("ABSENT")
        else:
            p(json.dumps(info, indent=2, default=str))

    p("\nLAYERS")
    if "layers" in f:
        layer_names = list(f["layers"].keys())
        p("Layer names:", layer_names)
        for name in layer_names:
            info = matrix_info(f, f"layers/{name}")
            p(json.dumps(info, indent=2, default=str))
    else:
        p("No layers group.")

    p("\nOBS AUDIT")
    if "obs" not in f:
        p("ERROR: obs absent")
    else:
        obs = f["obs"]
        p("obs shape/rows inferred from _index:",
          len(obs["_index"]) if "_index" in obs else "unknown")
        cols = [k for k in obs.keys() if k != "_index"]
        p("obs columns (n=%d):" % len(cols))
        p(cols)

        # Print value counts for likely biologically important metadata.
        wanted_exact = [
            "disease", "disease_ontology_term_id",
            "donor_id", "donor", "subject", "subject_id",
            "sex", "sex_ontology_term_id",
            "self_reported_ethnicity", "self_reported_ethnicity_ontology_term_id",
            "development_stage", "development_stage_ontology_term_id",
            "cell_type", "cell_type_ontology_term_id",
            "assay", "assay_ontology_term_id",
            "tissue", "tissue_ontology_term_id",
            "suspension_type", "is_primary_data",
            "dataset", "batch", "sample", "sample_id"
        ]
        # Also discover columns by keyword.
        keys = []
        for k in cols:
            kl = k.lower()
            if (k in wanted_exact or
                any(term in kl for term in [
                    "donor","subject","disease","diabet","condition",
                    "sex","ethnic","sample","batch","cell_type",
                    "development","assay","tissue"
                ])):
                keys.append(k)
        keys = list(dict.fromkeys(keys))

        p("\nSelected OBS value counts:")
        for k in keys:
            vals = read_vector(obs, k)
            if vals is None:
                p(f"\n[{k}] <not decoded by lightweight audit>")
                continue
            vc = value_counts(vals, top=60)
            p(f"\n[{k}] n={len(vals)} unique={len(np.unique(vals.astype(str)))}")
            for val, count in vc:
                p(f"  {val}: {count}")

    p("\nVAR AUDIT")
    if "var" not in f:
        p("ERROR: var absent")
    else:
        var = f["var"]
        nvar = len(var["_index"]) if "_index" in var else "unknown"
        p("var rows:", nvar)
        p("var columns:", [k for k in var.keys() if k != "_index"])
        if "_index" in var:
            idx = [dec(x) for x in var["_index"][:20]]
            p("first 20 var index values:", idx)

        for k in ["feature_name", "gene_symbols", "gene_symbol",
                  "feature_id", "gene_ids", "gene_id"]:
            if k in var:
                vals = read_vector(var, k)
                if vals is not None:
                    p(f"first 20 {k}:", [str(x) for x in vals[:20]])

    p("\nRAW VAR AUDIT")
    if "raw" in f and isinstance(f["raw"], h5py.Group) and "var" in f["raw"]:
        rv = f["raw/var"]
        p("raw/var rows:",
          len(rv["_index"]) if "_index" in rv else "unknown")
        p("raw/var columns:", [k for k in rv.keys() if k != "_index"])
        if "_index" in rv:
            p("first 20 raw/var index values:",
              [dec(x) for x in rv["_index"][:20]])
    else:
        p("raw/var absent")

p("\n" + "="*78)
p("AUDIT COMPLETE — NO EXPRESSION ANALYSIS PERFORMED")
p("="*78)

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")

print("\nSaved:", OUT)
