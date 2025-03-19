

fn(){

## 10X
input=(
`ls bigdata/10x/011725_NovaSeq-X_A/Tyler_Parse-10X_pool?/cellranger/outs/per_sample_outs/*/metrics_summary.csv`
)

for f in ${input[@]};do
    n=${f#*/per_sample_outs\/};n=${n%/*}
    cat $f | perl -ne 'if($_=~/([^,]+),\"?([\d,]+)(\s.+)?\"?$/){ print join("\t","'$n'_10x",$1,$2),"\n"; }' |\
    sed "s/Number of reads from cells called from this sample/number_of_reads/" |\
    sed "s/Total genes detected/number_of_genes/" |\
    sed "s/Mean reads per cell/mean_reads_per_cell/" |\
    sed "s/Median UMI counts per cell/median_umi_per_cell/" |\
    sed "s/Median genes per cell/median_genes_per_cell/" |\
    sed "s/Cells/number_of_cells/"  #|\
#    sed "s/UMIs per probe barcode/number_of_umi/" find matching id and extract next line
done  | grep -E "number|median|mean" | grep -v " "

## Parse
. src/util.sh
cat bigdata/parse/trailmaker/agg_samp_ana_summary.csv |\
python <( echo '
import sys
import pandas as pd
tt= pd.read_csv(sys.stdin)
y = tt.melt(id_vars=["statistic"], var_name="sample", value_name="Value")
y = y[["sample", "statistic", "Value"]]
y.to_csv(sys.stdout,index=False,header=False,sep="\t")
') | repl - <( cat bigdata/parse/manual/parse_seq_tymillerlab_nov2024.csv | grep primary_manual | tr "," "\t") 1 |\
sed "s/GRCh38-1-1-3c_//" |\
sed "s/tscp/umi/" |\
grep -E "median|mean|number" | sort -u  | perl -ne 'chomp; my@d=split/\t/,$_;$d[0]= (uc $d[0])."_parse"; print join("\t",@d),"\n";'
} 
fn | grep -v at50 | grep -v number_of_genes | perl -npe '$_=~s/\.\d+$//; $_ =~ s/(?<=\d)(?=(\d{3})+(?!\d))/,/g;' | python <( echo '
import sys
import pandas as pd
tt=pd.read_csv(sys.stdin,sep="\t",header=None)
tt.columns=["sample","metric","value"];
tt_wide = tt.pivot_table(index="metric", columns="sample", values="value", aggfunc="first")
tt_wide.to_csv(sys.stdout)
') > data/summary_metrics.csv
