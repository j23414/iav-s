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


// ==================== Equivalent to octoFLU ====================
process split_segments {
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

process uniq_merge {
  publishDir "results/", mode: 'symlink'
  input: path(clade_files)
  output: path("A0_metadata.tsv")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/nextstrain/dengue/new_ingest_uniqmerge/ingest/bin/uniq_merge.py
  chmod 755 uniq_merge.py
  python uniq_merge.py --cache H1_clade.txt --new H3_clade.txt --groupby_col barcode --outfile 01.txt
  python uniq_merge.py --cache 01.txt --new N1_clade.txt --groupby_col barcode --outfile 02.txt
  python uniq_merge.py --cache 02.txt --new N2_clade.txt --groupby_col barcode --outfile 03.txt
  python uniq_merge.py --cache 03.txt --new PB2_clade.txt --groupby_col barcode --outfile 04.txt
  python uniq_merge.py --cache 04.txt --new PB1_clade.txt --groupby_col barcode --outfile 05.txt
  python uniq_merge.py --cache 05.txt --new PA_clade.txt --groupby_col barcode --outfile 06.txt
  python uniq_merge.py --cache 06.txt --new M_clade.txt --groupby_col barcode --outfile 07.txt
  python uniq_merge.py --cache 07.txt --new NS_clade.txt --groupby_col barcode --outfile A0_metadata.tsv
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

process merge_motif{
  publishDir "results/", mode: 'symlink'
  input: tuple path(metadata),path(H3_motif)
  output: path("new_metadata.tsv")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/nextstrain/dengue/new_ingest_uniqmerge/ingest/bin/uniq_merge.py
  python uniq_merge.py --cache $metadata --new $H3_motif --groupby_col barcode --outfile new_metadata.tsv
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


workflow {
  print("Start")
  channel.fromPath(params.bvbrc_txt)
  | select_A0
  | select_GenBank_IDs
  | fetch_GenBank
  | genbank_to_fasta
  | split_segments
  | flatten
  | combine(genbank_to_fasta.out)
  | subset_fasta
  | Mafft
  | FastTree
  | treedist
  | get_clades
  | collect
  | uniq_merge
  | view

  /* pull h3 motif */
  h3_ch = uniq_merge.out
  | select_H3
  | view

  fetch_GenBank.out
  | genbank_to_protein
  | combine(select_H3.out)
  | select_H3_proteins
  | align_h3
  | get_h3_motif
  | view

  uniq_merge.out 
  | combine(get_h3_motif.out)
  | merge_motif
  | view
}