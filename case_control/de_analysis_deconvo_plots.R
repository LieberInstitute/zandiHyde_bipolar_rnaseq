library(jaffelab)
library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
library(GGally)
library(patchwork)
library(here)

load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint_deconvo.rda"), verbose = TRUE)
load(here("case_control","bipolarControl_deStats_byRegion_qSVAjoint.rda"), verbose = TRUE)

geneOut_deconvo2 <- geneOut_deconvo %>% 
  rownames_to_column("Gene") %>%
  as_tibble() %>%
  separate(Gene, into = c("term","Gene"), extra = "merge") %>%
  select(-deconvo_terms_Amyg, -deconvo_terms_sACC)

geneOut <- statOut %>%
  rownames_to_column("Gene") %>%
  as_tibble() %>%
  mutate(term = "no_deconvo") %>%
  filter(Gene %in% geneOut_deconvo2$Gene)

dim(geneOut)

## combine data
all_geneOut <- rbind(geneOut, geneOut_deconvo2)
dim(all_geneOut)

all_geneOut_sep <- map(c("Amyg", "sACC"), function(x){
  geneOut_region <- all_geneOut %>%
    select(Gene,term, ends_with(x)) %>%
    mutate(BrainRegion = x)
  
  colnames(geneOut_region) <- ss(colnames(geneOut_region),"_")
  return(geneOut_region)
})

all_geneOut_long <- do.call("rbind", all_geneOut_sep) %>%
  mutate(signif = cut(adj.P.Val, breaks = c(1, 0.05, 0.01, 0),
                      labels = c("<= 0.01", "<= 0.05", "> 0.05")))

all_geneOut_long %>% count(BrainRegion, term, signif)

## Compare t-stats
(signif_df <- all_geneOut_long  %>%
  select(Gene, signif, term, BrainRegion) %>%
  filter(signif != "> 0.05") %>%
  count(BrainRegion, signif, term) %>%
  pivot_wider(names_from = "term", values_from = "n"))

# BrainRegion signif    ilr no_deconvo  prop
# <chr>       <fct>   <int>      <int> <int>
#   1 Amyg        <= 0.01    25          5    33
# 2 Amyg        <= 0.05   233         99   291
# 3 sACC        <= 0.01    83        169    61
# 4 sACC        <= 0.05   296        497   257

#### ggpairs plots ####
t_wide <- all_geneOut_long %>% 
  select(Gene, BrainRegion, term, t) %>%
  pivot_wider(values_from = "t", names_from = "term")

t_ggpairs <- t_wide %>%
  group_by(BrainRegion) %>%
  group_map(~ggpairs(.x, c("no_deconvo", "prop","ilr"),
                     lower = list(continuous = wrap("smooth", alpha = 0.1, size=0.1, color = "blue"))))

names(t_ggpairs) <- t_wide %>%
  group_by(BrainRegion) %>%
  count() %>%
  pull(BrainRegion)

plot_fn <- here("deconvolution","plots","DE")

walk2(t_ggpairs, names(t_ggpairs), ~ggsave(.x + labs(title = .y),
                                          filename = paste0(plot_fn, "_ggpairs_tstats-",.y ,".png")))

## plot t-stats deconvo vs. no-deconvo 
t_stats <- all_geneOut_long %>%
  select(Gene, term, BrainRegion, t) %>%
  filter(term != "no_deconvo") %>%
  rename(deconvo.t = t) %>%
  left_join(all_geneOut_long %>%
              filter(term == "no_deconvo") %>%
              select(Gene, BrainRegion, t) %>%
              rename(no_deconvo.t = t) )
  
p_vals <- all_geneOut_long %>%
  select(Gene, term, BrainRegion, adj.P.Val) %>%
  filter(term != "no_deconvo") %>%
  rename(deconvo.adj.P.Val = adj.P.Val) %>%
  left_join(all_geneOut_long %>%
              filter(term == "no_deconvo") %>%
              select(Gene, BrainRegion, adj.P.Val) %>%
              rename(no_deconvo.adj.P.Val = adj.P.Val)) %>%
  mutate(Signif = case_when(no_deconvo.adj.P.Val < 0.05 & deconvo.adj.P.Val < 0.05 ~"sig_Both",
                                     no_deconvo.adj.P.Val < 0.05 ~ "sig_no-deconvo",
                                     deconvo.adj.P.Val < 0.05 ~ "sig_deconvo",
                                     TRUE ~ "None"))

t_limits <- all_geneOut_long %>%
  mutate(Neg = ifelse(t < 0, -1, 1)) %>%
  filter(signif == "<= 0.01") %>%
  group_by(BrainRegion, term, signif, Neg) %>%
  summarize(min_t = min(abs(t))) %>%
  mutate(min_t = min_t * Neg)

t_stats2 <- left_join(t_stats, p_vals)

signif_colors <- list(None = 'grey80', 
                      sig_deconvo = 'dark orange', 
                      `sig_no-deconvo` = 'skyblue3', 
                      sig_Both = 'purple')

t_lims_deconvo <- t_limits %>% filter(term != "no_deconvo")%>% ungroup()
t_lims_no_deconvo <- t_limits %>% filter(term == "no_deconvo")%>% ungroup() %>% select(-term) 

t_stat_scatter_both <- ggplot(t_stats2, aes(no_deconvo.t, deconvo.t, color = Signif))+
  geom_hline(data = t_lims_deconvo, aes(yintercept = min_t, linetype = signif))+
  geom_vline(data = t_lims_no_deconvo, aes(xintercept = min_t, linetype = signif))+
  geom_point(size = 0.5, alpha = 0.5) +
  facet_grid(term ~ BrainRegion) +
  scale_color_manual(values = signif_colors) +
  labs(x = "t-stat no deconvolution", y = "t-stat with deconvolution", 
       color = "FDR <= 0.05", linetype = "FDR <= 0.01")+
  theme_bw()

ggsave(t_stat_scatter_both, filename = paste0(plot_fn, "_t_stat_scatter_both.png"), width = 10 )

t_stat_scatter <- t_stats2 %>% filter(term == "prop") %>%
  ggplot(aes(no_deconvo.t, deconvo.t, color = Signif))+
  geom_hline(data = t_lims_deconvo %>% filter(term == "prop"), aes(yintercept = min_t, linetype = signif))+
  geom_vline(data = t_lims_no_deconvo, aes(xintercept = min_t, linetype = signif))+
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(~ BrainRegion, nrow = 1) +
  scale_color_manual(values = signif_colors) +
  labs(x = "t-stat no deconvolution", y = "t-stat with deconvolution", 
       color = "FDR <= 0.05", linetype = "FDR <= 0.01")+
  theme_bw()

ggsave(t_stat_scatter, filename = paste0(plot_fn, "_t_stat_scatter.pdf"), width = 10 )

t_stat_scatter_simple <- t_stats2 %>% filter(term == "prop") %>%
  ggplot(aes(no_deconvo.t, deconvo.t))+
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(~ BrainRegion, nrow = 1) +
  scale_color_manual(values = signif_colors) +
  labs(x = "t-stat no deconvolution", y = "t-stat with deconvolution")+
  theme_bw()

ggsave(t_stat_scatter_simple, filename = paste0(plot_fn, "_t_stat_scatter_simple.pdf"), width = 10 )
