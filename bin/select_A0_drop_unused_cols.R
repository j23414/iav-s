#! /usr/bin/env Rscript

library(tidyverse)

data <- read_delim("BVBRC_strain.txt", delim="\t")
names(data)
cdata <- data %>%
  filter(str_detect(strain, 'A0')) %>%
  select(-c(isolation_country, host_common_name, host_name, host_group, taxon_lineage_names, 
            public, species, genus, status, taxon_id, geographic_group, id, date_inserted, 
            date_modified, '_version_', segment_count, family))

write_delim(cdata, "BVBRC_metadata.tsv", delim="\t")
