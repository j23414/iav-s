#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.bvbrc_txt=false  // BVBRC_strain.txt

// === Import modules
include { select_A0;
          select_GenBank_IDs; } from './modules/bvbrc.nf'

include { split_Files;
          combine_Files; } from './modules/helper_functions.nf'

include { fetch_GenBank;
          genbank_to_fasta;
          genbank_to_protein; } from './modules/genbank.nf'

include { split_segments;
          subset_fasta;
          Mafft;
          FastTree;
          treedist;
          get_clades; 
          OctoFLU; } from './modules/octoflu.nf'

include { select_H3;
          select_H3_proteins;
          align_h3;
          get_h3_motif; } from './modules/h3motif.nf'

// === drafting processes
process cache_updatedate {
  publishDir "results/", mode: 'symlink'
  input: path(genbank_file)
  output: path("last_updated.txt")
  script:
  """
  #! /usr/bin/env bash
  cat $genbank_file \
    | grep "^LOCUS" \
    | awk 'OFS="\t" {print \$2,\$8}' \
    > last_updated.txt
  # e.g. "MZ892449        29-AUG-2021"
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

process uniq_merge {
  maxForks 8
  publishDir "results/", mode: 'symlink'
  input: path(clade_files)
  output: path("A0_metadata.tsv")
  script:
  """
  #! /usr/bin/env bash
  wget https://raw.githubusercontent.com/nextstrain/dengue/new_ingest_uniqmerge/ingest/bin/uniq_merge.py
  chmod 755 uniq_merge.py

  echo "barcode\tcollection_date\tH1_gb\tH3_gb\tN1_gb\tN2_gb\tPB2_gb\tPB1_gb\tPA_gb\tNP_gb\tM_gb\tNS_gb\tstrain\tus_clades\tH1\tH3\tN1\tN2\tPB2\tPB1\tPA\tNP\tM\tNS" > template.tsv
  echo "template_barcode\tcollection_date\tH1_gb\tH3_gb\tN1_gb\tN2_gb\tPB2_gb\tPB1_gb\tPA_gb\tNP_gb\tM_gb\tNS_gb\tstrain\tus_clades\tH1\tH3\tN1\tN2\tPB2\tPB1\tPA\tNP\tM\tNS" >> template.tsv
  
  python uniq_merge.py --cache template.tsv --new H1_clade.txt --groupby_col barcode --outfile 00.txt
  python uniq_merge.py --cache 00.txt --new H3_clade.txt --groupby_col barcode --outfile 01.txt
  python uniq_merge.py --cache 01.txt --new N1_clade.txt --groupby_col barcode --outfile 02.txt
  python uniq_merge.py --cache 02.txt --new N2_clade.txt --groupby_col barcode --outfile 03.txt
  rm 00.txt 01.txt 02.txt
  python uniq_merge.py --cache 03.txt --new PB2_clade.txt --groupby_col barcode --outfile 04.txt
  python uniq_merge.py --cache 04.txt --new PB1_clade.txt --groupby_col barcode --outfile 05.txt
  python uniq_merge.py --cache 05.txt --new PA_clade.txt --groupby_col barcode --outfile 06.txt
  rm 03.txt 04.txt 05.txt
  python uniq_merge.py --cache 06.txt --new M_clade.txt --groupby_col barcode --outfile 07.txt
  python uniq_merge.py --cache 07.txt --new NS_clade.txt --groupby_col barcode --outfile 08.txt
  cat 08.txt | grep -v "template_barcode" > A0_metadata.tsv
  rm 06.txt 07.txt 08.txt
  """
}

process merge_motif{
  maxForks 8
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

// === Main workflow
workflow {
  print("Start")
  channel.fromPath(params.bvbrc_txt)
  | select_A0
  | select_GenBank_IDs
  | split_Files // cache every 500 genbank entries
  | flatten
  | fetch_GenBank
  | collect
  | map{ n -> [n]}
  | combine(channel.of("A0.gb"))
  | combine_Files
  | genbank_to_fasta
  | OctoFLU
  // | split_segments
  // | flatten
  // | combine(genbank_to_fasta.out)
  // | subset_fasta
  // | Mafft
  // | FastTree
  // | treedist
  // | get_clades
  | collect
  | uniq_merge

  /* pull h3 motif */
  h3_ch = uniq_merge.out
  | select_H3

  combine_Files.out
  | genbank_to_protein
  | combine(select_H3.out)
  | select_H3_proteins
  | align_h3
  | get_h3_motif

  uniq_merge.out 
  | combine(get_h3_motif.out)
  | merge_motif
  | view

  /* thinking of caching */
  combine_Files.out
  | cache_updatedate
  | view
}