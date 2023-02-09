#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.bvbrc_txt=false  // BVBRC_strain.txt

process select_A0 {
  publishDir "results/", mode: 'symlink'
  input: path(bvbrc)
  output: path("bvbrc_metadata.tsv")
  script:
  """
  #! /usr/bin/env bash
  select_A0_drop_unused_cols.R
  """
}

process select_GenBank_IDs {
  publishDir "results/", mode: 'symlink'
  input: path(bvbrc)
  output: path("GenBank.ids")
  script:
  """
  #! /usr/bin/env bash
  cat $bvbrc \
    | tsv-select -H -f genbank_accessions \
    | grep -v "genbank_accessions" \
    | tr ';' '\n' \
    > GenBank.ids
  """
}

process select_New_Genbank_IDs {
  publishDir "results/", mode: 'symlink'
  input: tuple path(genbank_ids), path(cache_gb)
  output: path("new_GB.ids")
  script:
  """
  #! /usr/bin/env bash
  cat cache_gb \
    | grep "^LOCUS" \
    | awk -F'\t' '{print \$2}' \
    > new_GB.ids
  """
}

process fetch_GenBank {
  maxForks 1
  publishDir "results/", mode: 'symlink'
  input: path(genbank_ids)
  output: path("A0.gb")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/j23414/mini_nf/main/bin/batchFetchGB.sh
  chmod +x batchFetchGB.sh
  ./batchFetchGB.sh $genbank_ids > A0.gb
  """
}

process genbank_to_fasta {
  publishDir "results/", mode: 'symlink'
  input: path(genbank)
  output: path("A0.fasta")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/j23414/mini_nf/main/bin/procGenbank.pl
  chmod +x procGenbank.pl
  ./procGenbank.pl $genbank > A0.fasta
  """
}

// process run_octoFLU {
//   publishDir "results/", mode: 'symlink'
//   input: path(fasta)
//   output: path("octoflu/${fasta}_output")
//   script:
//   """
//   #! /usr/bin/env bash
//   git clone https://github.com/flu-crew/octoFLU.git
//   cd octoFLU
//   ./octoFLU.sh ../$fasta
//   """
// }

// process split_segments {
//   publishDir "results/", mode: 'symlink'
//   input: path(fasta)
//   output: path("*.ids")
//   script:
//   """
//   #! /usr/bin/env bash
//   # ===== Create your Blast Database
//   wget https://raw.githubusercontent.com/flu-crew/octoFLU/master/reference_data/reference.fa
// 
//   makeblastdb -in reference.fa -dbtype nucl
//   
//   # ===== Search your Blast Database
//   blastn -db reference.fa \
//     -query $fasta \
//     -num_alignments 1 \
//     -outfmt 6 \
//     -out blast_output.txt
//   
//   # ===== Split out query into 8 segments
//   cat blast_output.txt |awk -F'\t' '\$2~/\|H1\|/  {print \$1}' > H1.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|H3\|/  {print \$1}' > H3.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|N1\|/  {print \$1}' > N1.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|N2\|/  {print \$1}' > N2.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|PB2\|/ {print \$1}' > PB2.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|PB1\|/ {print \$1}' > PB1.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|PA\|/  {print \$1}' > PA.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|NP\|/  {print \$1}' > NP.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|M\|/   {print \$1}' > M.ids
//   cat blast_output.txt |awk -F'\t' '\$2~/\|NS\|/  {print \$1}' > NS.ids
//   """
// }

workflow {
  print("Start")
  channel.fromPath(params.bvbrc_txt)
  | select_A0
  | select_GenBank_IDs
  | fetch_GenBank
  | genbank_to_fasta
  | view
  //| split_segments
  //| view
}