
library(tidyverse)
library(here)
library(sessioninfo)
library(EnhancedVolcano)


load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)

statOut_long <- as.data.frame(statOut) %>% 
  select(starts_with("logFC"),starts_with("P.Val"),Type, Symbol) %>%
  rownames_to_column("mol") %>%
  pivot_longer(!c(mol,Type, Symbol)) %>% 
  separate(name,sep = "_", into = c("var","Region")) %>% 
  pivot_wider(names_from = "var", values_from = "value")

statOut_long$Type <- factor(statOut_long$Type, levels = c("Gene", "Exon", "Junction", "Transcript"))
head(statOut_long)

statOut_gene <- statOut_long %>% filter(Type == "Gene")

en_volcano <-  EnhancedVolcano(statOut_long,
                               lab = statOut_long$Symbol,
                                 x = 'logFC',
                                 y = 'P.Value') +
  facet_grid(Region~Type)

ggsave(en_volcano, filename = here("case_control","volcano_en.png"), width = 20, height = 10)

# sgejobs::job_single('de_results_compare_regions', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript de_results_compare_regions.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
