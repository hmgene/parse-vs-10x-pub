scanpy-marker(){
usage="$FUNCNAME <input.h5ad> [<marker.tsv]"
if [ $# -lt 1 ];then echo "$usage";return;fi 
echo '
import scanpy as sc
tt=sc.read_h5ad("'$1'")
mm = {
    "immune-cell": ["PTPRC"],
    "oligo": ["THY1", "PLP1", "APOD"],
    "pericyte": ["PDGFRB", "MCAM"],
    "endothelial": ["PECAM1", "ACTA2"],
    "cdc": ["HLA-DRA", "HLA-DRB1","HLA-DPB1", "AREG", "FCER1A","LAMP3",],
    "neutrophil": ["S100A8", "IFITM2", "FCGR3B"],
    "monocyte": ["LYZ", "VCAN", "FN1"],
    "b-cell": ["MS4A1"],
    "plasma": ["SDC1"],
    "mast-cell": ["MS4A2", "KIT"],
    "myeloid-cell": ["MPO","ITGAM","ITGAX","GPNMB", "TMEM119", "MRC1", "CD163", "IL1B", "CD83", "CCL3", "CD68","FUT4"],
    "t-cell": ["CD3E"],
    "glioma-cell": ["EGFR","PTPRZ1","CSRP2","AQP4","SOX2"]
}
sc.pl.dotplot(tt, mm1, groupby="cluster", standard_scale="var")

mm=pd.read_csv(m,sep="\t")
mm["genes"] = mm["genes"].str.rstrip(",").str.split(",")
mm1=mm.groupby("cell")["genes"].sum().to_dict()


# missing genes
mm.loc[~mm["gene"].isin(tt.var_names),:]
mm1=mm.loc[mm["gene"].isin(tt.var_names),:].groupby("cell-type")["gene"].apply(list).to_dict()
mm1 = mm.groupby("cell-type")["gene"].apply(list).to_dict()
sc.pl.dotplot(tt, mm1, groupby="cluster", standard_scale="var")
'

}
scanpy-10xdir2h5ad(){    
usage="$FUNCNAME <input_dir> <output.h5ad>"
if [ $# -lt 2 ];then echo "$usage";return;fi 
python <( echo '
import scanpy as sc
adata = sc.read_10x_mtx("'${1%/*}'"), var_names="gene_symbols", cache=True):
adata.write_h5ad("'${o%.h5ad}.h5ad'");
')
}

scanpy-prepro(){
usage="$FUNCNAME <input_adata.h5ad> <output_adata.h5ad>"
if [ $# -lt 2 ];then echo "$usage"; return; fi
i=${1:-"bigdata/output_combined/all-sample/DGE_filtered/anndata.h5ad"}
o=${2:-"bigdata/output_combined/all-sample/DGE_filtered/anndata_celltyped.h5ad"}
python <( echo '
i="'$i'";
o="'$o'";
import scanpy as sc
import numpy as np
tt0=sc.read_h5ad(i)
if "gene_name" in tt0.var :
    tt = tt0[:, ~tt0.var["gene_name"].duplicated()]
else:
    tt = tt0
tt.var.set_index("gene_name", inplace=True)
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
tt.write_h5ad(o)
')
}
