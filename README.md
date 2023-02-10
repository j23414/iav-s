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
[f6/790849] process > subset_fasta (6)       [100%] 10 of 10, cached: 10 ✔
[af/070bd5] process > Mafft (10)             [100%] 10 of 10, cached: 10 ✔
[15/c72d35] process > FastTree (10)          [100%] 10 of 10, cached: 10 ✔
[bd/700b0c] process > treedist (10)          [100%] 10 of 10, cached: 10 ✔
[23/8217dd] process > get_clades (10)        [100%] 10 of 10, cached: 10 ✔
[29/86ef37] process > uniq_merge             [100%] 1 of 1 ✔
[b0/8682ab] process > select_H3              [100%] 1 of 1 ✔
[45/c45eb9] process > genbank_to_protein (1) [100%] 1 of 1, cached: 1 ✔
[ef/897d82] process > select_H3_proteins (1) [100%] 1 of 1 ✔
[7c/4ab2f9] process > align_h3 (1)           [100%] 1 of 1 ✔
[f1/1f05a1] process > get_h3_motif (1)       [100%] 1 of 1 ✔
[bf/263290] process > merge_motif (1)        [100%] 1 of 1 ✔
/Users/jenchang/github/j23414/iav-s/work/bf/263290b08cb703cec00ac78fcf10bc/new_metadata.tsv

Completed at: 10-Feb-2023 08:45:52
Duration    : 10m 16s
CPU hours   : 1.2 (85.7% cached)
Succeeded   : 6
Cached      : 83
```

Thoughts on level of caching:

* Genbank level - most accurate, contains the entire GenBank entry, can get memory intensive, may be best genbank datasets up to 1GB
* Metadata level - smaller memory footprint, need this file anyway, may miss key information from GenBank
* ID level - smallest memory footprint, misses the "update date"
* Database - deduplicate sequences (shahash the sequence), faster access, maybe smaller memory foot print
* Maybe switch between levels of caching based on a flag `--cache-level 00_genbank` vs `--cache-level 01_metadata` with warning messages when datafiles get very large to indicate going up in cache level

Thoughts on merging:

* uniq merge takes a while for 10 segments across all strains
* could be parallized
* could only merge updated entries 
  * merge new data first 
  * subset cached metadata for update entries (by genbank id)
  * merge updated entries
  * concatinate with unchanged cached metadata
  * sort by collection date
  * ah, each segment's genbank has its own update date `H1_gb.update` and since each segment may have multiple genbanks...keep a separate "genbank\tlast_updated" file. 
