

import pandas as pd
import scanpy as sc

tt1=sc.read_h5ad("bigdata/10x/10x_merged.h5ad")
tt1.obs["platform"]="10x"
tt2=sc.read_h5ad("bigdata/parse/parse_merged.h5ad")
tt2.obs["platform"]="parse"

tt=sc.concat([tt1,tt2])
sc.pp.highly_variable_genes(tt, n_top_genes=4000, batch_key="group")
sc.pl.highly_variable_genes(tt)
sc.tl.pca(tt)
#sc.pl.pca_variance_ratio(tt, n_pcs=50, log=True)
sc.pp.neighbors(tt)
sc.tl.umap(tt)
sc.tl.leiden(tt,resolution=0.02, flavor="igraph", n_iterations=2)
sc.pl.umap( tt, color="sample", size=2)


