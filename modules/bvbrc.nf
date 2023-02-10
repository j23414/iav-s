#! /usr/bin/env nextflow

nextflow.enable.dsl=2

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