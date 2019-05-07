# blast2kraken
> generate kraken-style report from a blast results

kraken report是一种较直观地物种鉴定结果展示格式，基于LCA向上溯源，strain的reads_count结果会累加至species；

centrifuge-kreport 可将centrifuge结果生成kraken-style report; 在此脚本基础上，实现blast 结果生成kraken-style report脚本blast2kraken.pl的撰写；目的用于将nt库blast比对结果以更直观地方式展示物种鉴定结果，而不必再忧虑鉴定出一堆subgroup/strain物种;

Usage:
```
Usage: blast2kraken.pl -x <blast result>  -q <fasta>  -t <taxdb> OPTIONS  > <kraken-style.out>

blast2kraken.pl creates Kraken-style reports from blast out files.

Options:
    -x  Blast            (REQUIRED) Blast result
    -q  Fasta            (REQUIRED) Fasta input
    -t  TaxDB            (REQUIRED) Taxdb from taxonomy
    -min-ident  Score           Require a minimum identity score for reads alignment
    -min-length Score           Require a minimum lentgh for reads query
```
