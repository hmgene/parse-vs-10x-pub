import scanpy as sc
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

i="bigdata/parse/parse_merged.h5ad"
title="parse"
o="bigdata/parse/parse_merged"

tt=sc.read_h5ad(i)
r=pd.read_csv("data/meta-full.csv");
r1=r[["patient id","diagnosis"]]
r1.columns = ["group","diagnosis"] 
r1 = r1.drop_duplicates()
x = pd.merge(tt.obs, r1, on="group", how="left")
x.index = tt.obs.index 
tt.obs = x

m=pd.read_csv("data/cell_markers.tsv",sep="\t").set_index("cell_type")["markers"].str.split(",").to_dict()
for t,g in m.items():
    gg=[ i for i in g if i in tt.var.index ]
    if gg:
        tt.obs[t+"_score"] = tt[:,gg].X.mean(axis=1)

tt.obs["celltype"] = tt.obs[[t + "_score" for t in m]].idxmax(axis=1)

vp = sc.pl.stacked_violin(tt, [ i for i in tt.obs.keys() if "_score" in i],groupby="leiden",return_fig=True, swap_axes=True,figsize=(10,5),title="10x") 
vp.add_totals().style(ylim=(0,5))
vp.savefig(f"{o}_celltype_violin.png",dpi=300,bbox_inches="tight")


sc.pl.dotplot(tt,[ i for i in tt.obs.keys() if "_score" in i], groupby="leiden",swap_axes=True,show=False)
plt.savefig(f"{o}_celltype_dotplot.png",dpi=300, bbox_inches="tight")
plt.close()  # Close the figure

## umaps
sc.pl.umap(tt, color="leiden", size=2, show=False)
plt.savefig(f"{o}_celltype_umap.png", dpi=300, bbox_inches="tight")
plt.close()  # Close the figure

sc.pl.umap(tt, color="group", size=2, show=False)
plt.savefig(f"{o}_celltype_umap_group.png", dpi=300, bbox_inches="tight")
plt.close()  # Close the figure

sc.pl.umap(tt, color="diagnosis", size=2, show=False)
plt.savefig(f"{o}_celltype_umap_diagnosis.png", dpi=300, bbox_inches="tight")
plt.close()  # Close the figure

## end umaps
x = pd.crosstab(tt.obs["celltype"], tt.obs["group"])
p = x.div(x.sum(axis=0), axis=1)  # Normalize by column (Leiden cluster)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(8, 6))
sns.heatmap(x, annot=False, cmap='Blues', cbar_kws={'label': 'Score'}, xticklabels=x.columns, yticklabels=x.index)
#plt.xticks(rotation=45)  # Rotate x-axis labels
#plt.yticks(rotation=0)   # Rotate y-axis labels if needed
plt.title(title)
plt.tight_layout()  # Adjust layout to avoid label overlap
plt.savefig(f"{o}_celltype_leiden_heatmap.png")  # Saves to a file
plt.show()



import matplotlib.pyplot as plt
import pandas as pd

x = pd.crosstab(tt.obs["celltype"], tt.obs["leiden"])
fig, ax = plt.subplots(figsize=(10, 8))
row_totals = x.sum(axis=1)
ax.barh(x.index, row_totals, color='gray', alpha=0.7, height=0.8, label='Total', align='center')
ax.set_xlabel("Cells")
ax.set_ylabel("Cell Type")
ax.set_title(title)
plt.tight_layout()
plt.savefig(f"{o}_row_sum_bars.png", dpi=300, bbox_inches="tight")
#plt.show()












import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.imshow(x, cmap='coolwarm', interpolation='nearest')
plt.colorbar()  # To add color bar
plt.title( title)
plt.savefig(f"{o}_celltype_leiden_heatmap.png")  # Saves to a file
plt.close()  # Close the figure
#plt.show()

