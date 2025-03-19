

dotplot-marker(){
usage="$FUNCNAME <adata> <marker.csv> <outprefix>"
if [ $# -lt 2 ];then echo "$usage";return; fi
o=$3; mkdir -p ${o%/*};

python <( echo '
import scanpy as sc
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

i="'$1'";  #bigdata/10x/10x_merged.h5ad"
mf="'$2'"; #"data/cell_markers.tsv"
o="'$o'";  #bigdata/10x/10x_merged"

tt=sc.read_h5ad(i)
m=pd.read_csv(mf).set_index("cell_type")["markers"].str.split(",").to_dict()
for t,g in m.items():
    gg=[ i for i in g if i in tt.var.index ]
    if gg:
        tt.obs[t+"_score"] = tt[:,gg].X.mean(axis=1)

tt.obs["celltype"] = tt.obs[[t + "_score" for t in m]].idxmax(axis=1)

vp = sc.pl.stacked_violin(tt, [ i for i in tt.obs.keys() if "_score" in i],groupby="leiden",return_fig=True, swap_axes=True,figsize=(10,5),title="10x") 
vp.add_totals().style(ylim=(0,5))
vp.savefig(f"{o}_violin.png",dpi=300,bbox_inches="tight")
sc.pl.dotplot(tt,[ i for i in tt.obs.keys() if "_score" in i], groupby="leiden",swap_axes=True,show=False)
plt.savefig(f"{o}_dotplot.png",dpi=300, bbox_inches="tight")

')

}

