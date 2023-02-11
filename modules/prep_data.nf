#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process prep_upload {
  publishDir "results/upload", mode: 'symlink'
  input: tuple path(sequences), path(metadata)
  output: tuple path("*.tsv"), path("*.fasta")
  script:
  """
  #! /usr/bin/env bash
  
  for SEG in H1 H3 N1 N2 PB2 PB1 PA NP M NS ; do
    echo \$SEG
    cat $metadata \
    | tsv-filter --H --not-empty \${SEG}_gb \
    > \${SEG}_metadata.tsv

    cat \${SEG}_metadata.tsv \
    | tsv-select -H -f \${SEG}_gb \
    | grep -v "\${SEG}_gb" \
    > \${SEG}_ids.txt

    smof grep -f \${SEG}_ids.txt $sequences \
    > \${SEG}_sequences.fasta
  done

  # nextstrain remote upload remote_url metadata.tsv sequences.fasta
  """
}