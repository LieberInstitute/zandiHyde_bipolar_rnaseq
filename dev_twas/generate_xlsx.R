library(xlsx)
library(data.table)

load("xlsx_output.RData")

write.xlsx2(x = twas_z_amyg_threshold, file = "analysis/tables/BD_Amyg_sACC_FinalOutputTable.xlsx", sheetName = "Significant TWAS Z Scores in Amygdala", col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(x = twas_z_sacc_threshold, file = "analysis/tables/BD_Amyg_sACC_FinalOutputTable.xlsx", sheetName = "Significant TWAS Z Scores in sACC", col.names = TRUE, row.names = FALSE, append = TRUE)

write.xlsx2(x = twas_z_wide, file = "analysis/tables/BD_Amyg_sACC_FinalOutputTable.xlsx", sheetName = "TWAS Z Scatterplot with FDR and P-Values for Both Regions", col.names = TRUE, row.names = FALSE, append = TRUE)

write.xlsx2(x = merged_t, file = "analysis/tables/BD_Amyg_sACC_FinalOutputTable.xlsx", sheetName = "TWAS vs BD Differential Expression in Both Regions", col.names = TRUE, row.names = FALSE, append = TRUE)
