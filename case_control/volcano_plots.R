
library(tidyverse)
library(EnhancedVolcano)
library(here)
library(sessioninfo)

load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)

statOut_long <- as.data.frame(statOut) %>% 
  select(starts_with("logFC"),starts_with("adj.P.Val"),Type, Symbol) %>%
  rownames_to_column("mol") %>%
  pivot_longer(!c(mol,Type, Symbol)) %>% 
  separate(name,sep = "_", into = c("var","Region")) %>% 
  pivot_wider(names_from = "var", values_from = "value")

statOut_long$Type <- factor(statOut_long$Type, levels = c("Gene", "Exon", "Junction", "Transcript"))
head(statOut_long)

statOut_long %>% group_by(Region, Type) %>% 
  filter(abs(logFC) >1, adj.P.Val < 0.05) %>%
  summarize(n = n(),
            n_symb = paste(unique(Symbol), collapse = ", "))
# # A tibble: 6 x 4
# # Groups:   Region [2]
# Region Type         n n_symb              
# <chr>  <fct>    <int> <chr>               
#   1 Amyg   Gene         1 SLC5A5              
# 2 Amyg   Exon        14 IL1RL1, NR4A2, NR4A1
# 3 Amyg   Junction     1 NR4A1               
# 4 sACC   Gene         1 SERPINA3            
# 5 sACC   Exon         2 SLC11A1, HAMP       
# 6 sACC   Junction     2 HAMP, SLC11A1    

statOut_long %>%
  group_by(Region, Type) %>%
  summarise(max_lfc = max(logFC),
          min_lfc = min(logFC))

en_volcano <-  EnhancedVolcano(statOut_long,
                               lab = statOut_long$Symbol,
                                 x = 'logFC',
                                 y = 'adj.P.Val',
                               pCutoff = 0.05) +
  facet_grid(Region~Type) +
  xlim(-1.7, 1.7)+
  ylim(0, 5)+
  theme_bw(base_size = 10)

ggsave(en_volcano, filename = here("case_control","volcano_en.png"), width = 20, height = 10)
ggsave(en_volcano, filename = here("case_control","volcano_en.pdf"), width = 20, height = 10)

#### Simple Volcano ####

simple_volcano <- ggplot(statOut_long, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 0.05) + 
  facet_grid(Region~Type)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "lightskyblue3")+
  theme_bw(base_size = 10) +
  labs(x = "log2 Fold Change", y = "-log10 adj. P Value")

ggsave(simple_volcano, filename = here("case_control","volcano_simple.png"), width = 10, height = 5)
ggsave(simple_volcano, filename = here("case_control","volcano_simple.pdf"), width = 10, height = 5)

# sgejobs::job_single('volcano_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript volcano_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
