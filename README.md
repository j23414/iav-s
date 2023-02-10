# iav-s

glance at new api

Filter:

* host common name = Pig
* Isolation country = USA
* Keyword = A0

BV-BRC

* https://www.bv-brc.org/view/Taxonomy/11320#view_tab=strains_orthomyxoviridae&filter=and(eq(host_common_name,Pig),eq(isolation_country,%22USA%22),keyword(A0))

Download as Text `~/Downloads/BVBRC_strain.txt`

```
nextflow run main.nf --bvbrc_txt data/BVBRC_strain.txt

N E X T F L O W  ~  version 22.10.4
Launching `main.nf` [desperate_einstein] DSL2 - revision: c017b4b5a8
Start
executor >  local (1)
[a1/8746b3] process > select_A0 (1)          [100%] 1 of 1, cached: 1 ✔
[6a/105494] process > select_GenBank_IDs (1) [100%] 1 of 1, cached: 1 ✔
[ec/041846] process > fetch_GenBank (1)      [100%] 1 of 1, cached: 1 ✔
[e5/381521] process > genbank_to_fasta (1)   [100%] 1 of 1, cached: 1 ✔
[73/fda3ef] process > split_segments (1)     [100%] 1 of 1, cached: 1 ✔
[14/d38f8a] process > subset_fasta (9)       [100%] 10 of 10, cached: 10 ✔
[5e/32492b] process > Mafft (9)              [100%] 10 of 10, cached: 10 ✔
[ce/20537f] process > FastTree (10)          [100%] 10 of 10, cached: 10 ✔
[d6/386067] process > treedist (9)           [100%] 10 of 10, cached: 10 ✔
[2d/943644] process > get_clades (10)        [100%] 10 of 10, cached: 10 ✔
[5f/55e7ab] process > uniq_merge             [100%] 1 of 1 ✔
/Users/jchang3/github/j23414/iav-s/work/5f/55e7ab65fab2d001844c7131d48365/A0_metadata.tsv
Completed at: 09-Feb-2023 16:11:11
Duration    : 1m 16s
CPU hours   : 0.5 (95.9% cached)
Succeeded   : 1
Cached      : 55
```
