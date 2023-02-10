#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process split_Files {
  publishDir "results/", mode: 'symlink'
  input: path(largefile)
  output: path("A0_split*")
  script:
  """
  #! /usr/bin/env bash
  split -l 500 $largefile A0_split
  """
}

process combine_Files {
  publishDir "results/", mode: 'symlink'
  input: tuple path(list_of_files), val(filename)
  output: path("$filename")
  script:
  """
  #! /usr/bin/env bash
  cat $list_of_files > $filename
  """
}