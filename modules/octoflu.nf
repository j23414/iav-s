#! /usr/bin/env nextflow

nextflow.enable.dsl=2

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


// ==================== Equivalent to octoFLU ====================
process split_segments {
  maxForks 8
  publishDir "results/", mode: 'symlink'
  input: path(fasta)
  output: path("*.ids")
  script:
  """
  #! /usr/bin/env bash

  wget https://raw.githubusercontent.com/flu-crew/octoFLU/master/reference_data/reference.fa
  makeblastdb -in reference.fa -dbtype nucl
  
  blastn -db reference.fa \
    -query $fasta \
    -num_alignments 1 \
    -outfmt 6 \
    -out blast_output.txt
  
  cat blast_output.txt | awk -F'\t' '\$2~/\\|H1\\|/  {print \$1}' > H1.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|H3\\|/  {print \$1}' > H3.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|N1\\|/  {print \$1}' > N1.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|N2\\|/  {print \$1}' > N2.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|PB2\\|/ {print \$1}' > PB2.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|PB1\\|/ {print \$1}' > PB1.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|PA\\|/  {print \$1}' > PA.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|NP\\|/  {print \$1}' > NP.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|M\\|/   {print \$1}' > M.ids
  cat blast_output.txt | awk -F'\t' '\$2~/\\|NS\\|/  {print \$1}' > NS.ids
  """
}

process subset_fasta {
  maxForks 8
  publishDir "results/", mode: 'symlink'
  input: tuple path(ids), path(fasta)
  output: path("${ids.simpleName}.fasta")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/flu-crew/octoFLU/master/reference_data/reference.fa
  smof grep -Xf $ids $fasta > ${ids.simpleName}.fasta
  smof grep "|${ids.simpleName}|" reference.fa >> ${ids.simpleName}.fasta
  """
}

process Mafft {
  publishDir "results/", mode: 'symlink'
  input: path(fasta)
  output: path("${fasta.simpleName}_aln.fna")
  script:
  """
  #! /usr/bin/env bash
  mafft --auto $fasta > ${fasta.simpleName}_aln.fna
  """
}

process FastTree {
  publishDir "results/", mode: 'symlink'
  input: path(fasta)
  output: path("*.tre")
  script:
  """
  #! /usr/bin/env bash
  SEG=`echo "${fasta.simpleName}" | sed 's/_aln//g'`
  FastTree -nt $fasta > \${SEG}.tre
  """
}

process treedist {
  maxForks 8
  publishDir "results/", mode: 'symlink'
  input: path(tre)
  output: path("*.octoflu")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/flu-crew/octoFLU/master/treedist.py
  chmod 755 treedist.py

  SEG=`echo "${tre.simpleName}" | sed 's/_aln//g'`
  if [[ "\$SEG" == "H1" ]]; then
    COORDS="5,1,8"
  elif [[ "\$SEG" == "H3" ]]; then
    COORDS="5,1,8"
  else
    COORDS="5,1"
  fi

  python treedist.py -i $tre -c \$COORDS > \${SEG}.octoflu
  """
}

process get_clades {
  publishDir "results/", mode: 'symlink'
  input: path(octoflu)
  output: path("${octoflu.simpleName}_clade.txt")
  script:
  """
  #! /usr/bin/env bash
  if [[ "$octoflu.simpleName" == "H1" ]]; then
    echo "barcode\t${octoflu.simpleName}_gb\tstrain\tcollection_date\tgene\t${octoflu.simpleName}\tus_clades\texcess" > ${octoflu.simpleName}_clade.txt
  elif [[ "$octoflu.simpleName" == "H3" ]]; then
    echo "barcode\t${octoflu.simpleName}_gb\tstrain\tcollection_date\tgene\t${octoflu.simpleName}\tus_clades\texcess" > ${octoflu.simpleName}_clade.txt
  else
    echo "barcode\t${octoflu.simpleName}_gb\tstrain\tcollection_date\tgene\t${octoflu.simpleName}\texcess" > ${octoflu.simpleName}_clade.txt
  fi
  cat $octoflu \
    | awk -F'/' 'OFS="\t" {print \$4,\$0}'\
    | tr '|' '\t' \
    >> ${octoflu.simpleName}_clade.txt
  """
}

workflow OctoFLU {
  take:
    fasta_ch
  
  main:
    fasta_ch
      | split_segments
      | flatten
      | combine(fasta_ch)
      | subset_fasta
      | Mafft
      | FastTree
      | treedist
      | get_clades
  
  emit: get_clades.out
}