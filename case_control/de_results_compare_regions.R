
library(tidyverse)
library(here)
library(sessioninfo)

load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda"), verbose = TRUE)

head(statOut)
table(statOut$Type)

statOut <- as.data.frame(statOut) %>% 
  mutate(`FDR < 0.05` = case_when(adj.P.Val_sACC < 0.05 & adj.P.Val_Amyg < 0.05 ~"Both",
                                  adj.P.Val_sACC < 0.05 ~ "sACC",
                                  adj.P.Val_Amyg < 0.05 ~ "Amygdala",
                                  TRUE ~ "None"))

table(statOut$`FDR < 0.05`)
# Amygdala     Both     None     sACC 
#     253       67   759277     1768

statOut$Type <- factor(statOut$Type, levels = c("Gene", "Exon", "Junction", "Transcript"))

fdr_colors <- c(None = "gray", sACC = "skyblue3", Amygdala = "orange", Both = "purple")

t_stat_scatter <- ggplot(statOut, aes(x = t_Amyg, y = t_sACC, color = `FDR < 0.05`)) +
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(~Type) +
  labs(x = "t-stat Amygdala", y = "t-stat sACC") +
  scale_color_manual(values = fdr_colors) +
  theme_bw(base_size = 15)

ggsave(t_stat_scatter, filename = here("case_control","deStats_region_t_plots.pdf"))
ggsave(t_stat_scatter, filename = here("case_control","deStats_region_t_plots.png"))

# sgejobs::job_single('de_results_compare_regions', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript de_results_compare_regions.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
