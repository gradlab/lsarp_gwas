library(tidyverse)

arguments <- commandArgs(trailingOnly = T)

limits <- read_tsv(arguments[1], col_names = F)
p_threshold <- limits[[2,2]]

unitigs <- read_tsv(arguments[2], 
                    col_names = c("unitig", "af", "filter_p", "lrt_p", "beta", "beta_std_err", "h2", "annotation"))

unitigs <- unitigs %>% separate_wider_delim(annotation, delim = "\t", names = c("notes", "annotation"), too_few = "align_end")


unitigs <- unitigs %>% separate_wider_delim(annotation, 
                                            names = c("annotation1", "annotation2", "annotation3", "annotation4", "annotation5", "annotation6"), 
                                            delim = ',',
                                            too_few = "align_start") %>%
  pivot_longer(annotation1:annotation6, names_to = "annotation_number", values_to = "annotation", values_drop_na = T) %>%
  separate(annotation, c("reference", "location"), sep = ":") %>%
  separate(location, c("coordinates", "gene"), sep = ";") %>%
  separate(coordinates, c("start", "stop"), sep = "-") %>%
  mutate(start = as.numeric(start))

p <- ggplot(unitigs %>% mutate(significant = if_else(lrt_p < p_threshold, T, F)), 
       aes(x = start, y = -log10(lrt_p), color = significant)) +
  geom_point() + 
  scale_color_manual(values = c("#D3D3D3", arguments[3])) +
  theme_minimal() +
  xlab("Genomic position") +
  ylab("Significance (-log10(p))") +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  theme(legend.position = "none")

ggsave(arguments[4], p, width = 6, height = 4, units = "in", device="pdf")


