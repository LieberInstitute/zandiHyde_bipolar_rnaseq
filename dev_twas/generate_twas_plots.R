library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(plotly)
library(htmlwidgets)
library(sessioninfo)

# Sourcing Data/Inst. Vars. ####
load("rda/twas_exp_ranges.Rdata")

dir.create(file.path("analysis/plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("analysis/tables"), showWarnings = FALSE, recursive = TRUE)

# Filter N/A Z scores
twas_z <- twas_exp_fin %>% filter(!is.na(TWAS.Z))

twas_z_sacc <- twas_z[twas_z$region == "sacc", ]

twas_z_amyg <- twas_z[twas_z$region == "amygdala", ]

don <- list()

axisdf <- list()

don_key <- list()

p <- list()

intctv_plot <- list()

fin_plot <- list()

# don[[1]]$region
# don[[2]]$region

# Preprocessing Data ####
for (i in 1:2) {
    if (i == 1) {
        twas_var <- twas_z_amyg
    } else{
        twas_var <- twas_z_sacc
    }

    don[[i]] <-
        twas_var %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len = max(end)) %>%

        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(twas_var,
            .,
            by = c("CHR" = "CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, twas_mean_dist) %>%
        mutate(BPcum = twas_mean_dist + tot)


    axisdf[[i]] = don[[i]] %>% group_by(CHR) %>% summarise(center = (max(BPcum) + min(BPcum)) / 2)

    # Prepare text description for each SNP:
    don[[i]]$text <-
        paste0(
            "Gene Symbol: ",
            don[[i]]$genesymbol,
            "\nENSEMBL Gene ID: ",
            don[[i]]$geneid,
            "\nBrain Subregion: ",
            don[[i]]$region,
            "\nChromosome: ",
            don[[i]]$CHR,
            "\nStart Position: ",
            don[[i]]$start,
            "\nEnd Position: ",
            don[[i]]$end,
            "\nZ score: ",
            don[[i]]$TWAS.Z %>% round(2)
        )

    don_key[[i]] <-
        highlight_key(don[[i]], ~ genesymbol, group = "Gene Symbol")
}
# # Had to make sure these two were different
# test_thing_1 <- don_key[[1]]$data()
# test_thing_2 <- don_key[[2]]$data()
# md5(stri_paste(test_thing_1, collapse = ''))
# md5(stri_paste(test_thing_2, collapse = ''))
# i = 1

# TWAS Z Manhattan Plot ####
pdf(file = "analysis/plots/BIP_TWAS_ManhattanPlot.pdf")
# storing ggplot as an object3

sig <- qnorm(1 - 0.025 / table(twas_exp_fin$region))
for (i in 1:2) {
    # Bonferroni Correction
    sig_bonf <- sig[[i]]

    p[[i]] <-
        ggplot(don_key[[i]], aes(x = BPcum, y = TWAS.Z, text = text)) +

        ggtitle(paste0("Gene Windows of ", ifelse(i == 1, "Amygdala", "sACC") , " TWAS")) +
        # Show all points
        geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
        scale_color_manual(values = rep(c("#861657", "#D56AA0"), 22)) +
        geom_hline(
            yintercept = c(sig_bonf, -sig_bonf),
            color = "grey40",
            linetype = "dashed"
        ) +

        # custom X axis:
        scale_x_continuous(labels = axisdf[[i]]$CHR, breaks = axisdf[[i]]$center) +
        scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis

        # Custom the theme:
        theme_bw() +
        theme(
            legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        )

    print(p[[i]])
}
dev.off()

# Z scores threshold
twas_z_amyg_threshold <-
    rbind(twas_z_amyg[TWAS.Z > sig[[1]], ], twas_z_amyg[TWAS.Z < -sig[[1]], ])
twas_z_sacc_threshold <-
    rbind(twas_z_sacc[TWAS.Z > sig[[2]], ], twas_z_sacc[TWAS.Z < -sig[[2]], ])

twas_z_sig_tables <- list()


for (i in 1:2) {
    if (i == 1) {
        write.csv(twas_z_amyg_threshold, file = "analysis/tables/amygdala_twas_significant_genes_zscore.csv")
    } else {
        write.csv(twas_z_sacc_threshold, file = "analysis/tables/sacc_twas_significant_genes_zscore.csv")
    }
}

# Interactive TWAS Z Manhattan Plots ####
for (i in 1:2) {
    ##### Plotly
    intctv_plot[[i]] <- ggplotly(p[[i]], tooltip = "text")

    fin_plot[[i]] <- highlight(
        intctv_plot[[i]],
        on = "plotly_click",
        off = "plotly_doubleclick",
        color = "#60D394",
        selectize = TRUE
    )

    saveWidget(fin_plot[[i]],
        file.path(paste0(
            # "analysis/plots/",
            "BIP_TWAS_",
            ifelse(i == 1, "Amygdala", "sACC"),
            "_ManhattanPlotly.html"
        )))
}

# Issue #4 Plots ####
# https://github.com/LieberInstitute/zandiHyde_bipolar_rnaseq/issues/4

pdf('BIP_TWAS_ScatterPlots.pdf', useDingbats = FALSE, width = 10, height = 10)

# save.image(file = "ggplot_test.RData")
# load(file = "ggplot_test.RData")

# logical vector that indicates precense of gene in both subregions
twas_z[, in_both := uniqueN(region) == 2, by = c("start", "end")]

# create a column for each subregion where the values are z.score

# also include values where one region's z score is 0 for a gene

twas_z$fdr.z <- p.adjust(twas_z$TWAS.Z, 'fdr') # should be done by region not across all results, use this for ggplot
# Put into a table, csv, add it as a new sheet to the previous table you already made

twas_z <- select(twas_z, geneid, TWAS.Z, region) %>%
    as.data.table()

twas_z_wide <- dcast(twas_z, geneid ~ region, value.var = "TWAS.Z")

twas_z_wide$in_both <- ifelse(!is.na(twas_z_wide$amygdala & twas_z_wide$sacc), TRUE, FALSE)

# change this to reflect significant fdr < 0.05 z scores
twas_z_wide$FDR.5perc <- 'None'
twas_z_wide$FDR.5perc[twas_z_wide$amygdala < 0.05] <- 'amygdala'
twas_z_wide$FDR.5perc[twas_z_wide$sacc < 0.05] <- 'sACC'
twas_z_wide$FDR.5perc[twas_z_wide$amygdala < 0.05 & twas_z_wide$sacc < 0.05] <- 'Both'

twas_z_wide$FDR.5perc <- factor(twas_z_wide$FDR.5perc, levels = c('None', 'amygdala', 'sACC', 'Both'))

ggplot(twas_z_wide,
    aes(
        x = amygdala,
        y = sacc,
        color = FDR.5perc, # has four categories
        shape = in_both
     )) +
     geom_point() +
    coord_fixed() +
    theme_bw() +
    ggtitle('TWAS Z by brain region') +
    scale_color_manual(values = c('grey80', 'dark orange', 'skyblue3', 'purple')) # you can define names

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2020-11-13 17:04:12 EST"
# > proc.time()
#     user   system  elapsed
#    53.14     9.85 25140.28
# > options(width = 120)
# > session_info()
# - Session info -------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 4.0.3 (2020-10-10)
#  os       Windows 10 x64
#  system   x86_64, mingw32
#  ui       RStudio
#  language (EN)
#  collate  English_United States.1252
#  ctype    English_United States.1252
#  tz       America/New_York
#  date     2020-11-13
#
# - Packages -----------------------------------------------------------------------------------------------------------
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.0.2)
#  backports     1.1.10  2020-09-15 [1] CRAN (R 4.0.2)
#  callr         3.5.1   2020-10-13 [1] CRAN (R 4.0.3)
#  cli           2.1.0   2020-10-12 [1] CRAN (R 4.0.3)
#  colorspace    1.4-1   2019-03-18 [1] CRAN (R 4.0.2)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 4.0.2)
#  crosstalk     1.1.0.1 2020-03-13 [1] CRAN (R 4.0.3)
#  data.table  * 1.13.2  2020-10-19 [1] CRAN (R 4.0.3)
#  desc          1.2.0   2018-05-01 [1] CRAN (R 4.0.2)
#  devtools    * 2.3.2   2020-09-18 [1] CRAN (R 4.0.3)
#  digest        0.6.25  2020-02-23 [1] CRAN (R 4.0.2)
#  dplyr       * 1.0.2   2020-08-18 [1] CRAN (R 4.0.2)
#  ellipsis      0.3.1   2020-05-15 [1] CRAN (R 4.0.2)
#  fansi         0.4.1   2020-01-08 [1] CRAN (R 4.0.2)
#  farver        2.0.3   2020-01-16 [1] CRAN (R 4.0.2)
#  fastmap       1.0.1   2019-10-08 [1] CRAN (R 4.0.3)
#  fs            1.5.0   2020-07-31 [1] CRAN (R 4.0.2)
#  generics      0.1.0   2020-10-31 [1] CRAN (R 4.0.3)
#  ggplot2     * 3.3.2   2020-06-19 [1] CRAN (R 4.0.3)
#  ggrepel     * 0.8.2   2020-03-08 [1] CRAN (R 4.0.2)
#  glue          1.4.2   2020-08-27 [1] CRAN (R 4.0.2)
#  gtable        0.3.0   2019-03-25 [1] CRAN (R 4.0.2)
#  htmltools     0.5.0   2020-06-16 [1] CRAN (R 4.0.2)
#  htmlwidgets * 1.5.2   2020-10-03 [1] CRAN (R 4.0.3)
#  httpuv        1.5.4   2020-06-06 [1] CRAN (R 4.0.3)
#  httr          1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
#  jsonlite      1.7.1   2020-09-07 [1] CRAN (R 4.0.2)
#  labeling      0.4.2   2020-10-20 [1] CRAN (R 4.0.3)
#  later         1.1.0.1 2020-06-05 [1] CRAN (R 4.0.2)
#  lazyeval      0.2.2   2019-03-15 [1] CRAN (R 4.0.2)
#  lifecycle     0.2.0   2020-03-06 [1] CRAN (R 4.0.2)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 4.0.2)
#  memoise       1.1.0   2017-04-21 [1] CRAN (R 4.0.2)
#  mime          0.9     2020-02-04 [1] CRAN (R 4.0.0)
#  munsell       0.5.0   2018-06-12 [1] CRAN (R 4.0.2)
#  pillar        1.4.6   2020-07-10 [1] CRAN (R 4.0.2)
#  pkgbuild      1.1.0   2020-07-13 [1] CRAN (R 4.0.2)
#  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.0.2)
#  pkgload       1.1.0   2020-05-29 [1] CRAN (R 4.0.2)
#  plotly      * 4.9.2.1 2020-04-04 [1] CRAN (R 4.0.3)
#  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.0.2)
#  processx      3.4.4   2020-09-03 [1] CRAN (R 4.0.2)
#  promises      1.1.1   2020-06-09 [1] CRAN (R 4.0.2)
#  ps            1.3.4   2020-08-11 [1] CRAN (R 4.0.2)
#  purrr         0.3.4   2020-04-17 [1] CRAN (R 4.0.2)
#  R6            2.5.0   2020-10-28 [1] CRAN (R 4.0.3)
#  Rcpp          1.0.5   2020-07-06 [1] CRAN (R 4.0.2)
#  remotes       2.2.0   2020-07-21 [1] CRAN (R 4.0.2)
#  rlang         0.4.7   2020-07-09 [1] CRAN (R 4.0.2)
#  rprojroot     1.3-2   2018-01-03 [1] CRAN (R 4.0.2)
#  rstudioapi    0.11    2020-02-07 [1] CRAN (R 4.0.2)
#  scales        1.1.1   2020-05-11 [1] CRAN (R 4.0.2)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 4.0.2)
#  shiny         1.5.0   2020-06-23 [1] CRAN (R 4.0.3)
#  testthat    * 2.3.2   2020-03-02 [1] CRAN (R 4.0.2)
#  tibble        3.0.4   2020-10-12 [1] CRAN (R 4.0.3)
#  tidyr         1.1.2   2020-08-27 [1] CRAN (R 4.0.2)
#  tidyselect    1.1.0   2020-05-11 [1] CRAN (R 4.0.2)
#  usethis     * 1.6.3   2020-09-17 [1] CRAN (R 4.0.2)
#  utf8          1.1.4   2018-05-24 [1] CRAN (R 4.0.2)
#  vctrs         0.3.4   2020-08-29 [1] CRAN (R 4.0.2)
#  viridisLite   0.3.0   2018-02-01 [1] CRAN (R 4.0.2)
#  withr         2.3.0   2020-09-22 [1] CRAN (R 4.0.3)
#  xtable        1.8-4   2019-04-21 [1] CRAN (R 4.0.3)
#  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.0.0)
#
# [1] C:/Users/artas/Documents/R/win-library/4.0
# [2] C:/Program Files/R/R-4.0.3/library
