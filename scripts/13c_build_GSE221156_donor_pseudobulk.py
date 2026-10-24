#!/usr/bin/env python3
import os, sys, json, csv, hashlib
import numpy as np
try:
    import h5py
except ImportError:
    sys.exit("ERROR: h5py/numpy required. Activate project .venv first.")

H5 = "data/raw/GSE221156/GSE221156_Beta_processed_CELLxGENE.h5ad"
OUTDIR = "results/validation/GSE221156"
PROCDIR = "data/processed"
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(PROCDIR, exist_ok=True)

def decode_arr(x):
    a=np.asarray(x)
    if a.dtype.kind in ("S","O"):
        return np.array([v.decode("utf-8") if isinstance(v,(bytes,np.bytes_)) else str(v) for v in a], dtype=object)
    if a.dtype.kind=="U":
        return a.astype(object)
    return a

def read_anndata_col(group, name):
    obj=group[name]
    if isinstance(obj,h5py.Dataset):
        return decode_arr(obj[()])
    # categorical encoding
    keys=set(obj.keys())
    if "codes" in keys and "categories" in keys:
        codes=np.asarray(obj["codes"][()])
        cats=decode_arr(obj["categories"][()])
        out=np.empty(len(codes), dtype=object)
        for i,c in enumerate(codes):
            out[i]=None if c < 0 else cats[int(c)]
        return out
    # nullable/string arrays in newer encodings
    if "values" in keys:
        return decode_arr(obj["values"][()])
    raise ValueError(f"Unsupported obs encoding for {name}: keys={sorted(keys)}")

def read_index(group):
    idxname=group.attrs.get("_index","_index")
    if isinstance(idxname,bytes): idxname=idxname.decode()
    return decode_arr(group[idxname][()])

def fmt(x):
    if x is None: return ""
    if isinstance(x,float) and np.isnan(x): return ""
    return str(x)

print("="*78)
print("GSE221156 — DONOR METADATA AUDIT + RAW-COUNT PSEUDOBULK")
print("="*78)
print("Input:", H5)
if not os.path.exists(H5):
    sys.exit("ERROR: H5AD not found.")

with h5py.File(H5,"r") as f:
    obs=f["obs"]
    n_cells=len(read_index(obs))
    wanted=["donor_id","disease","sex","development_stage","BMI","BMI_Category",
            "HbA1c","Center","assay","self_reported_ethnicity"]
    cols={}
    for c in wanted:
        if c in obs:
            cols[c]=read_anndata_col(obs,c)
    for c,v in cols.items():
        if len(v)!=n_cells:
            sys.exit(f"ERROR: obs column length mismatch: {c}")

    donor=np.array([fmt(x) for x in cols["donor_id"]],dtype=object)
    disease=np.array([fmt(x) for x in cols["disease"]],dtype=object)
    donors=np.array(sorted(set(donor)),dtype=object)
    print(f"\nCells: {n_cells:,}")
    print(f"Donors: {len(donors)}")

    # donor-level metadata and consistency
    donor_rows=[]
    inconsistent=[]
    for d in donors:
        ix=np.where(donor==d)[0]
        row={"donor_id":d,"n_beta_cells":len(ix)}
        for c,v in cols.items():
            vals=[fmt(z) for z in v[ix] if fmt(z)!=""]
            uniq=sorted(set(vals))
            row[c] = ";".join(uniq)
            row[f"n_unique_{c}"]=len(uniq)
            if c in ["disease","sex","development_stage","Center","assay","self_reported_ethnicity"] and len(uniq)>1:
                # assay can legitimately vary, but record rather than fail
                if c!="assay":
                    inconsistent.append((d,c,uniq))
        donor_rows.append(row)

    # disease donor counts
    print("\nDONORS BY DISEASE")
    disease_donor_counts={}
    for row in donor_rows:
        dis=row["disease"]
        disease_donor_counts[dis]=disease_donor_counts.get(dis,0)+1
    for k in sorted(disease_donor_counts):
        print(f"  {k}: {disease_donor_counts[k]} donors")

    print("\nCELLS BY DISEASE")
    for k in sorted(set(disease)):
        print(f"  {k}: {int(np.sum(disease==k)):,} cells")

    if inconsistent:
        print("\nWARNING: donor-level metadata inconsistencies detected:")
        for x in inconsistent: print(" ",x)
    else:
        print("\nDonor-level disease/sex/age/center/ethnicity consistency: PASS")

    # save donor metadata
    meta_path=os.path.join(OUTDIR,"GSE221156_Donor_Metadata_Audit.csv")
    fields=[]
    for r in donor_rows:
        for k in r:
            if k not in fields: fields.append(k)
    with open(meta_path,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields)
        w.writeheader(); w.writerows(donor_rows)
    print("Saved:",meta_path)

    # feature symbols
    rv=f["raw/var"]
    if "feature_name" not in rv:
        sys.exit("ERROR: raw/var/feature_name absent.")
    genes=read_anndata_col(rv,"feature_name")
    genes=np.array([fmt(x) for x in genes],dtype=object)
    n_genes=len(genes)
    print(f"\nRaw genes: {n_genes:,}")
    print(f"Unique feature_name values: {len(set(genes)):,}")
    if len(set(genes)) != n_genes:
        print("WARNING: duplicate feature_name values exist; Ensembl IDs will also be retained.")

    ens=read_index(rv)
    ens=np.array([fmt(x) for x in ens],dtype=object)

    # raw sparse CSR
    X=f["raw/X"]
    shape=tuple(X.attrs["shape"])
    if shape!=(n_cells,n_genes):
        sys.exit(f"ERROR: raw/X shape {shape} != {(n_cells,n_genes)}")
    data=X["data"]; indices=X["indices"]; indptr=X["indptr"]
    print("raw/X shape:",shape)
    print("Stored nonzero values:",f"{len(data):,}")

    # donor index for each cell
    donor_to_i={d:i for i,d in enumerate(donors)}
    cell_donor_idx=np.array([donor_to_i[d] for d in donor],dtype=np.int16)

    # pseudobulk 48 x genes; float64 preserves source values exactly enough
    pb=np.zeros((len(donors),n_genes),dtype=np.float64)

    # Process cells in chunks, reading only corresponding CSR ranges.
    # Vectorized donor×gene bincount avoids a Python loop over ~337M nonzeros.
    chunk_cells=500
    print("\nAggregating raw counts by donor (chunked CSR; no dense cell×gene matrix)...")
    flat_size=len(donors)*n_genes
    for start in range(0,n_cells,chunk_cells):
        end=min(start+chunk_cells,n_cells)
        ptr=np.asarray(indptr[start:end+1],dtype=np.int64)
        a=int(ptr[0]); b=int(ptr[-1])
        vals=np.asarray(data[a:b],dtype=np.float64)
        idx=np.asarray(indices[a:b],dtype=np.int64)
        local_ptr=ptr-a
        nnz_per_cell=np.diff(local_ptr)
        donor_for_nz=np.repeat(cell_donor_idx[start:end].astype(np.int64),nnz_per_cell)
        flat_index=donor_for_nz*n_genes+idx
        chunk_sum=np.bincount(flat_index,weights=vals,minlength=flat_size)
        pb += chunk_sum.reshape(len(donors),n_genes)
        if start==0 or end==n_cells or end%10000==0:
            print(f"  processed {end:,}/{n_cells:,} cells")

    # integer-like audit of pseudobulk
    nz=pb[pb!=0]
    frac_integer=np.mean(np.isclose(nz,np.rint(nz),atol=1e-6)) if nz.size else np.nan
    print("\nPseudobulk nonzero integer-like fraction:",f"{frac_integer:.8f}")
    print("Total raw counts aggregated:",f"{pb.sum():.0f}")

    # Save TSV (genes rows, donors columns) for R/edgeR
    pb_path=os.path.join(PROCDIR,"GSE221156_Beta_DonorPseudobulk_RawCounts.tsv")
    print("Writing:",pb_path)
    with open(pb_path,"w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n")
        w.writerow(["ensembl_id","gene_symbol"]+donors.tolist())
        for g in range(n_genes):
            w.writerow([ens[g],genes[g]]+[str(int(round(x))) for x in pb[:,g]])

    # compact candidate audit directly from pseudobulk
    candidate_genes=[
        "SEMA7A","EDIL3","FLRT3","PGF","SPP1",
        "RPH3A","PPIA","FGF13","FGF12","SST","CNTN2","LGALS8","RSPO2",
        "ITGB1","ITGAV","UNC5B","UNC5D","NRP1","CD44","ITGB5","ITGA1","PLXNC1",
        "BSG","FGFR1","SSTR1","SSTR2","SSTR3","CNTN1","NRCAM","ITGA3","LGR4","RNF43"
    ]
    symbol_to_idx={}
    for i,g in enumerate(genes):
        symbol_to_idx.setdefault(g,[]).append(i)

    cand_path=os.path.join(OUTDIR,"GSE221156_CandidateGene_DonorPseudobulk_Audit.csv")
    cand_rows=[]
    lib=pb.sum(axis=1)
    for g in candidate_genes:
        ids=symbol_to_idx.get(g,[])
        if not ids:
            cand_rows.append({"gene":g,"present":False})
            continue
        counts=pb[:,ids].sum(axis=1)
        cpm=np.divide(counts,lib,out=np.zeros_like(counts),where=lib>0)*1e6
        for di,d in enumerate(donors):
            cand_rows.append({
                "gene":g,"present":True,"donor_id":d,
                "disease":next(r["disease"] for r in donor_rows if r["donor_id"]==d),
                "n_beta_cells":next(r["n_beta_cells"] for r in donor_rows if r["donor_id"]==d),
                "raw_count":int(round(counts[di])),
                "CPM":float(cpm[di])
            })
    cfields=["gene","present","donor_id","disease","n_beta_cells","raw_count","CPM"]
    with open(cand_path,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cfields)
        w.writeheader(); w.writerows(cand_rows)
    print("Saved:",cand_path)

# hashes
for fp in [meta_path,pb_path,cand_path]:
    h=hashlib.md5()
    with open(fp,"rb") as fh:
        for block in iter(lambda:fh.read(1024*1024),b""): h.update(block)
    print("MD5",os.path.basename(fp),h.hexdigest())

print("\nPRIMARY COMPARISON PRESPECIFIED:")
print("  T2D vs normal at DONOR level; prediabetes excluded from primary DE.")
print("  Prediabetes retained for secondary trajectory/direction analysis.")
print("  No cell-level pseudoreplication.")
print("\n"+"="*78)
print("SCRIPT 13c COMPLETE")
print("="*78)
