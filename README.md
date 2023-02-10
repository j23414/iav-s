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
[5e/32492b] process > Mafft (10)             [100%] 10 of 10, cached: 10 ✔
[83/e13da5] process > FastTree (9)           [100%] 10 of 10, cached: 10 ✔
[15/afcaa6] process > treedist (10)          [100%] 10 of 10, cached: 10 ✔
[5b/d9c978] process > get_clades (10)        [100%] 10 of 10, cached: 10 ✔
[7f/43dc28] process > uniq_merge             [100%] 1 of 1, cached: 1 ✔
[af/b85566] process > select_H3              [100%] 1 of 1, cached: 1 ✔
[78/18b8af] process > genbank_to_protein (1) [100%] 1 of 1, cached: 1 ✔
[c9/30cbf4] process > select_H3_proteins (1) [100%] 1 of 1, cached: 1 ✔
[e3/103a9a] process > align_h3 (1)           [100%] 1 of 1, cached: 1 ✔
[64/52ad18] process > get_h3_motif (1)       [100%] 1 of 1, cached: 1 ✔
[10/99e7f4] process > merge_motif (1)        [100%] 1 of 1 ✔
/Users/jchang3/github/j23414/iav-s/work/7f/43dc28efbaf4357185d6add5adb130/A0_metadata.tsv
/Users/jchang3/github/j23414/iav-s/work/af/b8556653b17ea62229ac2b70ba00ff/H3.gb_ids
/Users/jchang3/github/j23414/iav-s/work/64/52ad1853b7bfd133e19a7f714197c5/H3_motif.txt
/Users/jchang3/github/j23414/iav-s/work/10/99e7f44d5f4af0e76c9682b0c54906/new_metadata.tsv
Completed at: 09-Feb-2023 16:11:11
Duration    : 1m 16s
CPU hours   : 0.5 (95.9% cached)
Succeeded   : 1
Cached      : 55
```
