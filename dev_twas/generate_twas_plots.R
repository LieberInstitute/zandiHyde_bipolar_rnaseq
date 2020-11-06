library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(plotly)
library(htmlwidgets)
library(openssl)
library(stringi)
library(sessioninfo)

load("rda/twas_exp_ranges.Rdata")

# Filter N/A Z scores
twas_z <- twas_exp_fin %>% filter(!is.na(TWAS.Z))

twas_z_sacc <- twas_z[twas_z$region == "sacc",]

twas_z_amyg <- twas_z[twas_z$region == "amygdala",]

don <- list()

axisdf <- list()

don_key <- list()

p <- list()

intctv_plot <- list()

fin_plot <- list()

# don[[1]]$region
# don[[2]]$region

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
pdf(file = "BIP_TWAS_ManhattanPlot.pdf")
# storing ggplot as an object3

# Bonferroni Correction?
# sig <- 0.05 / nrow(twas_exp_fin)
# i <- 1
for (i in 1:2) {
    p[[i]] <-
        ggplot(don_key[[i]], aes(x = BPcum, y = TWAS.Z, text = text)) +

        ggtitle(paste0("Gene Windows of ", ifelse(i == 1, "Amygdala", "sACC") , " TWAS")) +
        # Show all points
        geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
        scale_color_manual(values = rep(c("#861657", "#D56AA0"), 22)) +
        # geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
        # geom_hline(yintercept = log10(sig), color = "grey40", linetype = "dashed") +
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
        paste0(
            "BIP_TWAS_",
            ifelse(i == 1, "Amygdala", "sACC"),
            "_ManhattanPlotly.html"
        ))
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info

# > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2020-11-06 12:37:52 EST"
# > proc.time()
#    user  system elapsed
#   2.798   6.101 115.679
# > options(width = 120)
# > session_info()
# ??? Session info ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  setting  value
#  version  R version 4.0.2 Patched (2020-06-24 r78746)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-11-06
#
# ??? Packages ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  package     * version date       lib source
#  askpass       1.1     2019-01-13 [2] CRAN (R 4.0.0)
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 4.0.0)
#  cli           2.1.0   2020-10-12 [2] CRAN (R 4.0.2)
#  colorout    * 1.2-2   2020-09-24 [1] Github (jalvesaq/colorout@726d681)
#  colorspace    1.4-1   2019-03-18 [2] CRAN (R 4.0.0)
#  crayon        1.3.4   2017-09-16 [2] CRAN (R 4.0.0)
#  data.table  * 1.13.2  2020-10-19 [2] CRAN (R 4.0.2)
#  digest        0.6.26  2020-10-17 [2] CRAN (R 4.0.2)
#  dplyr       * 1.0.2   2020-08-18 [2] CRAN (R 4.0.2)
#  ellipsis      0.3.1   2020-05-15 [2] CRAN (R 4.0.0)
#  fansi         0.4.1   2020-01-08 [2] CRAN (R 4.0.0)
#  generics      0.0.2   2018-11-29 [2] CRAN (R 4.0.0)
#  ggplot2     * 3.3.2   2020-06-19 [1] CRAN (R 4.0.2)
#  ggrepel     * 0.8.2   2020-03-08 [2] CRAN (R 4.0.0)
#  glue          1.4.2   2020-08-27 [2] CRAN (R 4.0.2)
#  gtable        0.3.0   2019-03-25 [2] CRAN (R 4.0.0)
#  htmltools     0.5.0   2020-06-16 [2] CRAN (R 4.0.2)
#  htmlwidgets * 1.5.2   2020-10-03 [2] CRAN (R 4.0.2)
#  httr          1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
#  jsonlite      1.7.0   2020-06-25 [1] CRAN (R 4.0.2)
#  lazyeval      0.2.2   2019-03-15 [2] CRAN (R 4.0.0)
#  lifecycle     0.2.0   2020-03-06 [2] CRAN (R 4.0.0)
#  magrittr      1.5     2014-11-22 [2] CRAN (R 4.0.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 4.0.0)
#  openssl     * 1.4.2   2020-06-27 [1] CRAN (R 4.0.2)
#  pillar        1.4.6   2020-07-10 [2] CRAN (R 4.0.2)
#  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.0.0)
#  plotly      * 4.9.2.1 2020-04-04 [1] CRAN (R 4.0.2)
#  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.0.0)
#  R6            2.4.1   2019-11-12 [2] CRAN (R 4.0.0)
#  Rcpp          1.0.5   2020-07-06 [1] CRAN (R 4.0.2)
#  rlang         0.4.8   2020-10-08 [1] CRAN (R 4.0.2)
#  scales        1.1.1   2020-05-11 [2] CRAN (R 4.0.0)
#  sessioninfo * 1.1.1   2018-11-05 [2] CRAN (R 4.0.0)
#  stringi     * 1.5.3   2020-09-09 [2] CRAN (R 4.0.2)
#  tibble        3.0.4   2020-10-12 [2] CRAN (R 4.0.2)
#  tidyr         1.1.2   2020-08-27 [2] CRAN (R 4.0.2)
#  tidyselect    1.1.0   2020-05-11 [2] CRAN (R 4.0.0)
#  vctrs         0.3.4   2020-08-29 [2] CRAN (R 4.0.2)
#  viridisLite   0.3.0   2018-02-01 [2] CRAN (R 4.0.0)
#  withr         2.3.0   2020-09-22 [2] CRAN (R 4.0.2)
#
# [1] /users/aseyedia/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
