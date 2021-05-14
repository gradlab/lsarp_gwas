library(tidyverse)

arguments <- commandArgs(trailingOnly = T)

# read in annotated, significant unitigs and rename columns
annotated <- read_tsv(arguments[1], col_names = F)
colnames(annotated) <- c("unitig", "af", "filter_p", "lrt_p", "beta", "beta_std_err", "h2", "notes", "annotation")
annotated <- annotated %>% select(-filter_p, -notes)

# filter out unitigs that map to multiple locations
annotated_uniquemap <- annotated %>% 
    filter(!str_detect(annotation, ",")) %>%
    arrange(lrt_p)

# separate reference genome into a new column
annotated_uniquemap <- annotated_uniquemap %>%
    separate(annotation, c("reference", "location"), sep = ":")

# count the reference genome chosen for the top 10% most significant unitigs

annotated_uniquemap %>% 
    top_frac(-0.1, lrt_p) %>%
    count(reference) %>%
    top_n(1, n) %>%
    pull(reference)
