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
[e9/82ad03] process > select_A0 (1)          [100%] 1 of 1, cached: 1 ✔
[90/e15900] process > select_GenBank_IDs (1) [100%] 1 of 1, cached: 1 ✔
[a5/f1511d] process > split_Files (1)        [100%] 1 of 1, cached: 1 ✔
[13/b033ac] process > fetch_GenBank (26)     [100%] 26 of 26, cached: 26 ✔
[6d/b6c3e9] process > combine_Files (1)      [100%] 1 of 1, cached: 1 ✔
[c3/c71518] process > genbank_to_fasta (1)   [100%] 1 of 1, cached: 1 ✔
[60/1f203a] process > split_segments (1)     [100%] 1 of 1, cached: 1 ✔
[61/c01ede] process > subset_fasta (10)      [100%] 10 of 10, cached: 10 ✔
[74/20f53e] process > Mafft (10)             [100%] 10 of 10, cached: 10 ✔
[de/e16d73] process > FastTree (10)          [100%] 10 of 10, cached: 10 ✔
[5c/b37035] process > treedist (10)          [100%] 10 of 10, cached: 10 ✔
[54/7fc575] process > get_clades (10)        [100%] 10 of 10, cached: 10 ✔
[b2/6cbe55] process > uniq_merge             [  0%] 0 of 1
[-        ] process > select_H3              -
[45/c45eb9] process > genbank_to_protein (1) [100%] 1 of 1, cached: 1 ✔
[-        ] process > select_H3_proteins     -
[-        ] process > align_h3               -
[-        ] process > get_h3_motif           -
[-        ] process > merge_motif            -
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
