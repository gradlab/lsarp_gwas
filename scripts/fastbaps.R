library(fastbaps)
library(phytools)
library(readr)

snp_fasta <- import_fasta_sparse_nt("data/gubbins/201008-tm__gubbins/core_alignment.filtered_polymorphic_sites.fasta", 
    prior = "baps")

gubbins_tree <- read.newick("data/gubbins/201008-tm__gubbins/core_alignment.final_tree.tre")
gubbins_tree_rooted <- midpoint.root(gubbins_tree)

best_partition <- best_baps_partition(snp_fasta, gubbins_tree_rooted)

best_partition_df <- data.frame(id = gubbins_tree_rooted$tip.label, fastbaps = best_partition)

write_tsv(best_partition_df, "data/population_structure/fastbaps.txt")
