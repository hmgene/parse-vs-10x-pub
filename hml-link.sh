#data/10x/pool1_UH14_web_summary.html
x="https://raw.githack.com/hmgene/parse-vs-10x-pub/main/"
for f in data/*/*_*_*.html;do
    a=${f#*data/};a=${a%/*}
    b=${f%\_summary*};b=${b%_*};b=${b##*/}
    o="[$b]($x$f)"
    echo $o
done
