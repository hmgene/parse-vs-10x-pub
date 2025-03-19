import scanpy as sc
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

o="bigdata/parse/merged"
tt=sc.read_h5ad("bigdata/parse/scanpy/anndata.h5ad")
tt.obs.merge(pd.read_csv("bigdata/parse/manual/parse_seq_tymillerlab_nov2024.csv"))
tmp=tt.obs.merge(pd.read_csv("bigdata/parse/manual/parse_seq_tymillerlab_nov2024.csv"),on="sample",how="left") 
tmp.index = tt.obs.index
tt.obs=tmp

tt = tt[tt.obs["condition"] == "primary_manual", :]
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
sc.pp.highly_variable_genes(tt, n_top_genes=4000, batch_key="group")
sc.pl.highly_variable_genes(tt)
sc.tl.pca(tt)
sc.pl.pca_variance_ratio(tt, n_pcs=50, log=True)
sc.pp.neighbors(tt)
sc.tl.umap(tt)
sc.tl.leiden(tt,resolution=0.02, flavor="igraph", n_iterations=2)
sc.pl.umap( tt, color="leiden", size=2)

m=pd.read_csv("data/cell_markers.tsv",sep="\t").set_index("cell_type")["markers"].str.split(",").to_dict()

for t,g in m.items():
    m[t]=[ i for i in g if i in tt.var.index ]


tt.write_h5ad(f"{o}.h5ad");
p=sc.pl.umap(tt, color="leiden", size=2, show=False)
plt.savefig(f"{o}_umap_leiden.png", dpi=300, bbox_inches="tight")

p = sc.pl.dotplot(tt,m, groupby="leiden",show=False)
plt.savefig(f"{o}_dotplot_leiden.png", dpi=300, bbox_inches="tight")


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

