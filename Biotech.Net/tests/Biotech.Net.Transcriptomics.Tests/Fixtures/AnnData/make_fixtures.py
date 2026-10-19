"""Generates the AnnData/h5ad fixtures and scanpy reference values used by AnnDataTests and AtlasScaleTests.

Run with a Python environment that has anndata, scanpy, h5py, scikit-misc and scikit-learn:

    python make_fixtures.py [output_dir]

Everything is seeded and the package versions are written to versions.txt. The pbmc3k subset is downloaded
by scanpy on first run (10x Genomics "3k PBMCs from a healthy donor", a public dataset).
"""
import gzip
import importlib.metadata as md
import os
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

out = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(20261019)


def save_csv(name, df, fmt="%.17g"):
    with gzip.open(os.path.join(out, name), "wt", newline="") as f:
        df.to_csv(f, index=True, float_format=fmt)


# ---------------------------------------------------------------- 1. kitchen sink, dense X
n, m = 24, 15
X = rng.poisson(2.0, size=(n, m)).astype(np.float32)
obs = pd.DataFrame(index=pd.Index([f"cell{i}" for i in range(n)], name="barcode"))
obs["cluster"] = pd.Categorical(rng.choice(["T", "B", "NK"], n), categories=["B", "NK", "T"])
obs["stage"] = pd.Categorical(rng.choice(["early", "mid", "late"], n), categories=["early", "mid", "late"], ordered=True)
obs["batch_num"] = pd.Categorical(rng.choice([1, 2], n))  # numeric categories
obs["n_counts"] = X.sum(1).astype(np.float64)
obs["n_genes"] = (X > 0).sum(1).astype(np.int32)
obs["score"] = rng.normal(size=n).astype(np.float32)
obs["big"] = np.arange(n, dtype=np.int64) * 10_000_000_000
obs["is_doublet"] = rng.random(n) < 0.2
obs["donor"] = [f"donor_{i % 4}" for i in range(n)]
ni = pd.array(rng.integers(0, 50, n), dtype="Int64")
ni[[1, 5]] = pd.NA
obs["nullable_int"] = ni
nb = pd.array(rng.random(n) < 0.5, dtype="boolean")
nb[[0, 7]] = pd.NA
obs["nullable_bool"] = nb
obs["with_missing_cat"] = pd.Categorical(["a", None, "b"] * 8)
var = pd.DataFrame(index=[f"GENE{j}" for j in range(m)])
var["gene_symbol"] = [("MT-" if j < 2 else "") + f"SYM{j}" for j in range(m)]
var["highly_variable"] = rng.random(m) < 0.5
var["means"] = X.mean(0).astype(np.float64)
a = ad.AnnData(X=X, obs=obs, var=var)
a.obsm["X_pca"] = rng.normal(size=(n, 5))
a.obsm["X_umap"] = rng.normal(size=(n, 2)).astype(np.float32)
a.varm["PCs"] = rng.normal(size=(m, 5))
a.obsp["distances"] = sp.random(n, n, density=0.2, format="csr", random_state=1, dtype=np.float64)
a.varp["corr"] = sp.random(m, m, density=0.3, format="csr", random_state=2, dtype=np.float32)
a.layers["counts"] = sp.csr_matrix(X.astype(np.int32))
a.layers["dense_layer"] = X * 2
a.uns["title"] = "kitchen sink"
a.uns["version"] = 3
a.uns["pi"] = 3.25
a.uns["flag"] = True
a.uns["ints"] = np.array([1, 2, 3], dtype=np.int64)
a.uns["names"] = np.array(["x", "y", "zed"], dtype=object)
a.uns["matrix"] = np.arange(6, dtype=np.float64).reshape(2, 3)
a.uns["cluster_colors"] = ["#1f77b4", "#ff7f0e", "#2ca02c"]
a.uns["neighbors"] = {"params": {"n_neighbors": 15, "method": "umap", "metric": "euclidean"}, "connectivities_key": "connectivities"}
a.uns["nested"] = {"level1": {"level2": {"value": 1.5, "text": "deep"}}}
a.raw = a.copy()
a.write_h5ad(os.path.join(out, "kitchen_sink.h5ad"))

# ---------------------------------------------------------------- 2. CSR X, gzip-compressed and chunked
Xs = sp.random(40, 30, density=0.25, format="csr", random_state=3, dtype=np.float32)
Xs.data = np.round(Xs.data * 10).astype(np.float32) + 1
b = ad.AnnData(X=Xs, obs=pd.DataFrame(index=[f"c{i}" for i in range(40)]), var=pd.DataFrame(index=[f"g{j}" for j in range(30)]))
b.obs["group"] = pd.Categorical(["g1", "g2"] * 20)
b.write_h5ad(os.path.join(out, "csr_gzip.h5ad"), compression="gzip")

# ---------------------------------------------------------------- 3. CSC X, float64
c = ad.AnnData(X=sp.csc_matrix(Xs.astype(np.float64)), obs=b.obs.copy(), var=b.var.copy())
c.write_h5ad(os.path.join(out, "csc.h5ad"))

# Exact values for the reader tests, dumped independently of any HDF5 reader.
import json  # noqa: E402

with open(os.path.join(out, "expected_values.json"), "w") as f:
    json.dump({
        "kitchen_X": X.astype(np.float64).ravel().tolist(),
        "kitchen_X_pca": a.obsm["X_pca"].ravel().tolist(),
        "kitchen_X_umap": a.obsm["X_umap"].astype(np.float64).ravel().tolist(),
        "kitchen_score": a.obs["score"].astype(np.float64).tolist(),
        "kitchen_big": a.obs["big"].tolist(),
        "kitchen_cluster": [str(v) for v in a.obs["cluster"]],
        "kitchen_is_doublet": [bool(v) for v in a.obs["is_doublet"]],
        "kitchen_nullable_int": [None if pd.isna(v) else int(v) for v in a.obs["nullable_int"]],
        "kitchen_nullable_bool": [None if pd.isna(v) else bool(v) for v in a.obs["nullable_bool"]],
        "kitchen_distances": a.obsp["distances"].toarray().ravel().tolist(),
        "kitchen_PCs": a.varm["PCs"].ravel().tolist(),
        "csr_dense": Xs.toarray().astype(np.float64).ravel().tolist(),
    }, f)

# ---------------------------------------------------------------- 4. legacy layout (pre-0.7 h5sparse X, dataframe 0.1.0 with __categories)
with h5py.File(os.path.join(out, "legacy_v01.h5ad"), "w") as f:
    small = sp.csr_matrix(Xs[:6, :5])
    g = f.create_group("X")
    g.attrs["h5sparse_format"] = "csr"
    g.attrs["h5sparse_shape"] = np.array(small.shape)
    g["data"] = small.data
    g["indices"] = small.indices
    g["indptr"] = small.indptr
    vlen = h5py.special_dtype(vlen=str)
    for name, idx, cols in (("obs", [f"c{i}" for i in range(6)], {"label": ["u", "v", "u", "w", "v", "u"]}),
                            ("var", [f"g{j}" for j in range(5)], {})):
        d = f.create_group(name)
        d.attrs["encoding-type"] = "dataframe"
        d.attrs["encoding-version"] = "0.1.0"
        d.attrs["_index"] = "index"
        d.attrs["column-order"] = np.array(list(cols), dtype=vlen)
        d.create_dataset("index", data=np.array(idx, dtype=object), dtype=vlen)
        if cols:
            cats = d.create_group("__categories")
            for k, vals in cols.items():
                cat = pd.Categorical(vals)
                cd = cats.create_dataset(k, data=np.array(cat.categories, dtype=object), dtype=vlen)
                cd.attrs["ordered"] = False
                codes = d.create_dataset(k, data=cat.codes.astype(np.int8))
                codes.attrs["categories"] = cd.ref
legacy = ad.read_h5ad(os.path.join(out, "legacy_v01.h5ad"))
assert list(legacy.obs["label"]) == ["u", "v", "u", "w", "v", "u"], legacy.obs

# ---------------------------------------------------------------- 5. pbmc3k subset + scanpy references
full = sc.datasets.pbmc3k()
full.var_names_make_unique()
sc.pp.filter_genes(full, min_cells=3)
cells = np.sort(rng.choice(full.n_obs, 300, replace=False))
# Genes detected in at least 5 of the chosen cells: rarer genes share identical (mean, variance) pairs and make
# scikit-misc's loess (seurat_v3) numerically singular on a subset this small.
detected = np.asarray((full.X[cells] > 0).sum(0)).ravel()
mt = np.where(full.var_names.str.startswith("MT-") & (detected >= 5))[0]
others = np.setdiff1d(np.where(detected >= 5)[0], mt)
genes = np.sort(np.concatenate([mt, rng.choice(others, 1000 - len(mt), replace=False)]))
p = full[cells, genes].copy()
p.X = sp.csr_matrix(p.X, dtype=np.float32)
p.obs = pd.DataFrame(index=p.obs_names)
p.var = pd.DataFrame({"gene_ids": full.var["gene_ids"].values[genes]}, index=p.var_names)
p.write_h5ad(os.path.join(out, "pbmc3k_subset.h5ad"), compression="gzip")

q = p.copy()
q.var["mt"] = q.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(q, qc_vars=["mt"], inplace=True, log1p=True)
save_csv("pbmc_qc_obs.csv.gz", q.obs)
save_csv("pbmc_qc_var.csv.gz", q.var.drop(columns=["gene_ids"]))

hv3 = sc.pp.highly_variable_genes(p, flavor="seurat_v3", n_top_genes=250, inplace=False)
save_csv("pbmc_hvg_seurat_v3.csv.gz", hv3)

lg = p.copy()
sc.pp.normalize_total(lg, target_sum=1e4)
sc.pp.log1p(lg)
med = p.copy()
sc.pp.normalize_total(med)  # median target
lgx = lg.X.tocsr()
save_csv("pbmc_lognorm_sums.csv.gz", pd.DataFrame({
    "row_sum": np.asarray(lgx.sum(1, dtype=np.float64)).ravel(),
    "median_norm_row_sum": np.asarray(med.X.sum(1, dtype=np.float64)).ravel()}, index=lg.obs_names))
save_csv("pbmc_lognorm_colsums.csv.gz", pd.DataFrame({"col_sum": np.asarray(lgx.sum(0, dtype=np.float64)).ravel()}, index=lg.var_names))
coo = lgx.tocoo()
pick = rng.choice(coo.nnz, 200, replace=False)
with open(os.path.join(out, "pbmc_lognorm_entries.csv"), "w", newline="") as f:
    f.write("row,col,value\n")
    for k in pick:
        f.write(f"{coo.row[k]},{coo.col[k]},{float(coo.data[k]):.9g}\n")

hvs = sc.pp.highly_variable_genes(lg, flavor="seurat", n_top_genes=250, inplace=False)
save_csv("pbmc_hvg_seurat.csv.gz", hvs)

hvg_idx = np.where(hv3["highly_variable"].values)[0]
sub = lg[:, hvg_idx].copy()
sc.pp.pca(sub, n_comps=20, zero_center=True, svd_solver="arpack", random_state=0)
dense = sub.X.toarray().astype(np.float64)
from sklearn.decomposition import PCA  # noqa: E402

skp = PCA(n_components=20, svd_solver="full").fit(dense)
save_csv("pbmc_pca_variance.csv.gz", pd.DataFrame({
    "scanpy_variance": sub.uns["pca"]["variance"], "scanpy_ratio": sub.uns["pca"]["variance_ratio"],
    "sklearn_variance": skp.explained_variance_, "sklearn_ratio": skp.explained_variance_ratio_}))
save_csv("pbmc_pca_loadings.csv.gz", pd.DataFrame(skp.components_[:10].T, index=sub.var_names), "%.10g")
save_csv("pbmc_pca_scores.csv.gz", pd.DataFrame(skp.transform(dense)[:, :10], index=sub.obs_names), "%.10g")

with open(os.path.join(out, "versions.txt"), "w", newline="\n") as f:
    for pkg in ("anndata", "scanpy", "h5py", "numpy", "scipy", "pandas", "scikit-learn", "scikit-misc"):
        f.write(f"{pkg} {md.version(pkg)}\n")
print("done")
