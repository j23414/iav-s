#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process select_H3 {
  publishDir "results/", mode: 'symlink'
  input: path(metadata)
  output: path("H3.gb_ids")
  script:
  """
  #! /usr/bin/env bash
  cat $metadata \
    | tsv-select -H -f H3_gb \
    | grep -v "^\$" \
    | grep -v "H3_gb" \
    | tr ',' '\n' \
    > H3.gb_ids
  """
}

process select_H3_proteins {
  publishDir "results/", mode: 'symlink'
  input: tuple path(protein_fasta), path(H3_ids)
  output: path("H3.faa")
  script:
  """
  #! /usr/bin/env bash
  smof grep -f $H3_ids $protein_fasta > H3.faa
  """
}

process align_h3 {
  publishDir "results/", mode: 'symlink'
  input: path(H3_faa)
  output: path("H3_aln.faa")
  script:
  """
  #! /usr/bin/env bash
  mafft --auto $H3_faa > H3_aln.faa
  """
}

process get_h3_motif {
  publishDir "results/", mode: 'symlink'
  input: path(H3_aln_faa)
  output: path("H3_motif.txt")
  script:
  """
  #! /usr/bin/env bash
  # A regular expression that should capture the entire HA1 region
  H1_HA_regex="DT[LI]C.*QSR"
  H3_HA_regex="QKL.*QTR"
  MIN_LENGTH=100
  bounds=\$(smof grep -qP --gff  "\$H1_HA_regex|\$H3_HA_regex" $H3_aln_faa | cut -f4,5 | sort | uniq -c | sort -rg | head -1 | sed 's/ *[0-9]* *//')
  smof subseq -b \$bounds $H3_aln_faa | smof clean -u > ${H3_aln_faa.simpleName}_HA1.faa

  # Bolton et al. 2019 "Antigenic evolution of H3N2 influenza A viruses in swine in the United States from 2012 to 2016"
  echo "barcode\tH3_motif" > H3_motif.txt
  motif.pl 145,155,156,158,159,189 ${H3_aln_faa.simpleName}_HA1.faa \
    | awk -F'/' 'OFS="\t" {print \$4,\$0}' \
    | awk -F'\t' 'OFS="\t" {print \$1,\$2}' >> H3_motif.txt
  """
}