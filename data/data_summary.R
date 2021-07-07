library(SummarizedExperiment)
library(tidyverse)

## load data
load("zandiHypde_bipolar_rseGene_n511.rda")
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Controls","Cases"),
                     levels = c("Cases", "Controls"))

rse_gene$BrainRegion <- factor(rse_gene$BrainRegion, levels = c("sACC", "Amygdala"))

pd <- as.data.frame(colData(rse_gene))
colnames(pd)

msd_sum <- function(value){
  m <- mean(value)
  s <- sd(value)

  return(paste0(sprintf("%.2f",m)," (",sprintf("%.2f",s),")"))
}

data_sum <- pd %>% group_by(BrainRegion, Dx) %>%
  summarise(n = n(),
            n_F = sum(Sex == "F"),
            `N Female (%)` = paste0(n_F, " (",sprintf("%.2f",n_F*100/n),"%)"),
            `Mean Age (SD)` = msd_sum(AgeDeath),
            `Mean PMI (SD)` = msd_sum(PMI),
            `Mean pH (SD)` = msd_sum(pH),
            `Mean RIN (SD)` = msd_sum(RIN)) %>%
  mutate(r = paste(BrainRegion, Dx)) %>%
  column_to_rownames("r") %>%
  select(-n_F, -BrainRegion, -Dx)

(data_sum2 <- t(data_sum))
# sACC Cases      sACC Controls   Amygdala Cases  Amygdala Controls
# n             "126"           "142"           "121"           "122"
# N Female (%)  "51 (40.48%)"   "24 (16.90%)"   "45 (37.19%)"   "26 (21.31%)"
# Mean Age (SD) "42.59 (12.90)" "50.47 (15.67)" "42.95 (13.24)" "52.30 (16.48)"
# Mean PMI (SD) "28.98 (13.42)" "28.93 (12.07)" "28.62 (13.84)" "29.95 (12.29)"
# Mean pH (SD)  "3.86 (3.10)"   "4.41 (3.02)"   "4.03 (3.07)"   "4.09 (3.09)"
# Mean RIN (SD) "7.90 (0.87)"   "7.78 (0.75)"   "7.55 (0.90)"   "7.29 (0.76)"

write.csv(data_sum2, "supptable_new1.csv")
