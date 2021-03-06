# blast2kraken
> generate kraken-style report from a blast results

kraken report是一种较直观地物种鉴定结果展示格式，基于 taxonomy-tree LCA溯源，strain的reads_count结果会累加至species level, species reads_count结果会累加至genus level，...

centrifuge-kreport 可将centrifuge结果生成kraken-style report; 在此脚本基础上，实现blast 结果生成kraken-style report脚本blast2kraken.pl的撰写；目的用于将nt库blast比对结果以更直观地方式展示物种鉴定结果，而不必再忧虑鉴定出一堆subgroup/strain物种;

> Usage:
```
Usage: blast2kraken.pl -x <blast result>  -q <fasta>  -t <taxdb> OPTIONS  > <kraken-style.out>
       blast2kraken.pl -x <blast result> -t <taxdb> OPTIONS > <kraken-style.out>

blast2kraken.pl creates Kraken-style reports from blast out files.

Options:
    -x  Blast            (REQUIRED) Blast result
    -t  TaxDB            (REQUIRED) Taxdb from taxonomy
    -q  Fasta		 Fasta input, to count all sequencing reads
    -min-ident  Score           Require a minimum identity score for reads alignment
    -min-length Score           Require a minimum lentgh for reads query
```
注：此脚本适用于blast -outfmt '6 qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname'格式的比对结果；

若blast格式不符，可更改脚本中对应项的位置；
```
my ($readID,$qlen,$seqID,$pident,$taxID)= @arr[0,1,2,5,15] ;
```
> Example结果鉴定物种展示
```
$less example.blast2kraken.report|awk -F '\t' '$1>1 && $4=="S"{sub(/^[ \t]+|[ ]+$/,"",$NF);print $1"\t"$2"\t"$4"\t"$5"\t"$NF}'|head -5
 23.98  29981   S       1280    Staphylococcus aureus
  4.00  4996    S       1639    Listeria monocytogenes
  3.98  4978    S       1423    Bacillus subtilis
  3.96  4946    S       1613    Lactobacillus fermentum
  3.95  4936    S       1351    Enterococcus faecalis
```
