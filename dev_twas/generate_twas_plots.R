library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(plotly)
library(htmlwidgets)
library(openssl)
library(stringi)

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
for (i in 1:2) {
    p[[i]] <-
        ggplot(don_key[[i]], aes(x = BPcum, y = TWAS.Z, text = text)) +

        ggtitle(paste0("Gene Windows of ", ifelse(i == 1, "Amygdala", "sACC") , " TWAS")) +
        # Show all points
        geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
        scale_color_manual(values = rep(c("#861657", "#D56AA0"), 22)) +

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
