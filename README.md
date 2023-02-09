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

N E X T F L O W  ~  version 21.10.6
Launching `main.nf` [distracted_einstein] - revision: 78693406d6
Start
executor >  local (1)
[24/6c3bdb] process > select_A0 (1) [100%] 1 of 1 ✔
/Users/jenchang/github/j23414/iav-s/work/24/6c3bdbf05457ee1b034b09179cfd27/bvbrc_metadata.tsv
```
