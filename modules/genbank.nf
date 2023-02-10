#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process fetch_GenBank {
  maxForks 1
  publishDir "results/", mode: 'symlink'
  input: path(genbank_ids)
  output: path("${genbank_ids.simpleName}.gb")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/j23414/mini_nf/main/bin/batchFetchGB.sh
  chmod +x batchFetchGB.sh
  ./batchFetchGB.sh $genbank_ids > ${genbank_ids.simpleName}.gb
  """
}

process genbank_to_fasta {
  maxForks 8
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

process genbank_to_protein {
  publishDir "results/", mode: 'symlink'
  input: path(genbank)
  output: path("${genbank.simpleName}.faa")
  script:
  """
  #! /usr/bin/env bash
  genbank2protein.pl $genbank > ${genbank.simpleName}.faa
  """
}
