import scanpy as sc
import pandas as pd 

tt1=sc.read_h5ad("bigdata/10x/011725_NovaSeq-X_A/Tyler_Parse-10X_pool1/scanpy/Tyler_Parse-10X_pool1.h5ad")
tt2=sc.read_h5ad("bigdata/10x/011725_NovaSeq-X_A/Tyler_Parse-10X_pool2/scanpy/Tyler_Parse-10X_pool2.h5ad")
m1=pd.read_csv("bigdata/10x/011725_NovaSeq-X_A/Tyler_Parse-10X_pool1/scanpy/meta.csv")
m2=pd.read_csv("bigdata/10x/011725_NovaSeq-X_A/Tyler_Parse-10X_pool2/scanpy/meta.csv")

tt=tt1.concatenate(tt2)
tt.var["mt"] = tt.var_names.str.startswith("MT-")
tt.obs["pct_counts_mt"] = np.sum(tt[:, tt.var["mt"]].X, axis=1).A1 / np.sum(tt.X, axis=1).A1 * 100
tt.var["ribo"] = tt.var_names.str.startswith(("RPS", "RPL"))
tt.var["hb"] = tt.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics( tt, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
sc.pp.filter_cells(tt, min_genes=200)
sc.pp.filter_genes(tt, min_cells=3)
tt= tt[tt.obs["pct_counts_mt"] < 25, :]
tt.layers["counts"] = tt.X.copy()
sc.pp.normalize_total(tt)
sc.pp.log1p(tt)
sc.pp.highly_variable_genes(tt, n_top_genes=4000, batch_key="sample")
sc.pl.highly_variable_genes(tt)
sc.tl.pca(tt)
sc.pl.pca_variance_ratio(tt, n_pcs=50, log=True)
sc.pp.neighbors(tt)
sc.tl.umap(tt)
sc.tl.leiden(tt,resolution=0.02, flavor="igraph", n_iterations=2)
sc.pl.umap( tt, color="leiden", size=2)
sc.pl.umap( tt, color="sample", size=2)


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
df = tt.obs.copy()
co_occurrence = pd.crosstab(df['sample'], df['leiden'])
x=np.log10(co_occurrence+1)
plt.figure(figsize=(10, 10))
#sns.heatmap(np.log10(co_occurrence+1), annot=True, cmap="Blues", fmt="d", linewidths=0.5)
sns.clustermap(x, method="ward", metric="euclidean", cmap="coolwarm", annot=True)
plt.xlabel("Leiden Cluster")
plt.ylabel("Sample")
plt.title("Co-occurrence Heatmap: Leiden vs Sample")
plt.show()

