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
nextflow run j23414/iav-s -r main --bvbrc_txt data/BVBRC_strain.txt -resume

N E X T F L O W  ~  version 22.10.4
Launching `main.nf` [desperate_einstein] DSL2 - revision: c017b4b5a8
Start
executor >  local (1)
[e9/82ad03] process > select_A0 (1)              [100%] 1 of 1, cached: 1 ✔
[90/e15900] process > select_GenBank_IDs (1)     [100%] 1 of 1, cached: 1 ✔
[a5/f1511d] process > split_Files (1)            [100%] 1 of 1, cached: 1 ✔
[13/b033ac] process > fetch_GenBank (26)         [100%] 26 of 26, cached: 26 ✔
[6d/b6c3e9] process > combine_Files (1)          [100%] 1 of 1, cached: 1 ✔
[c3/c71518] process > genbank_to_fasta (1)       [100%] 1 of 1, cached: 1 ✔
[23/01d543] process > OctoFLU:split_segments (1) [100%] 1 of 1 ✔
[bb/c2d6b3] process > OctoFLU:subset_fasta (6)   [100%] 10 of 10 ✔
[dc/05741f] process > OctoFLU:Mafft (9)          [100%] 10 of 10 ✔
[14/38011c] process > OctoFLU:FastTree (8)       [100%] 10 of 10 ✔
[21/914d07] process > OctoFLU:treedist (10)      [100%] 10 of 10 ✔
[4e/7d15ae] process > OctoFLU:get_clades (10)    [100%] 10 of 10 ✔
[1b/ac29a6] process > uniq_merge                 [100%] 1 of 1 ✔
[e0/1c8f36] process > select_H3                  [100%] 1 of 1 ✔
[45/c45eb9] process > genbank_to_protein (1)     [100%] 1 of 1, cached: 1 ✔
[3e/e47730] process > select_H3_proteins (1)     [100%] 1 of 1 ✔
[2f/be63d6] process > align_h3 (1)               [100%] 1 of 1 ✔
[ca/42b4e9] process > get_h3_motif (1)           [100%] 1 of 1 ✔
[ad/d5b37c] process > merge_motif (1)            [100%] 1 of 1 ✔
/Users/jenchang/github/j23414/iav-s/work/ad/d5b37c7e17b8f4425917cd60f78e42/new_metadata.tsv

Completed at: 10-Feb-2023 09:52:13
Duration    : 18m 7s
CPU hours   : 0.8 (14.9% cached)
Succeeded   : 57
Cached      : 32
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

Thoughts on api limits:

* If pulling data from a database, set `maxForks 1` up to `maxForks 3`
* If pulling scripts from github, may not be necessary but could use `maxForks 8` to be safe
