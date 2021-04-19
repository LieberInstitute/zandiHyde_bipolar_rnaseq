
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

statOut$Type <- factor(statOut$Type, levels = c("Gene", "Exon", "Junction", "Transcript"))
statOut$`FDR < 0.05` <- factor(statOut$`FDR < 0.05`, levels = c("Both", "Amygdala", "sACC", "None"))

table(statOut$`FDR < 0.05`)
# Both Amygdala     sACC     None 
# 67      253     1768   759277 

## spearman cor
region_cor <- statOut %>% 
  group_by(Type) %>%
  summarize(cor = cor(t_Amyg, t_sACC, method = "spearman")) %>%
  mutate(cor_anno = paste0("rho==", format(round(cor, 2), nsmall = 2)))

fdr_colors <- c(None = "gray", sACC = "skyblue3", Amygdala = "orange", Both = "purple")

t_stat_scatter <- ggplot(statOut, aes(x = t_Amyg, y = t_sACC)) +
  geom_point(aes(color = `FDR < 0.05`), size = 0.5, alpha = 0.5) +
  facet_wrap(~Type) +
  labs(x = "t-stat Amygdala", y = "t-stat sACC") +
  scale_color_manual(values = fdr_colors) +
  geom_text(data = region_cor, aes(x = -4.25, y = 6, label = cor_anno), parse = TRUE) +
  theme_bw(base_size = 15) +
  NULL

ggsave(t_stat_scatter, filename = here("case_control","deStats_region_t_plots.png"))
ggsave(t_stat_scatter, filename = here("case_control","deStats_region_t_plots.pdf"))

# sgejobs::job_single('de_results_compare_regions', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript de_results_compare_regions.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
