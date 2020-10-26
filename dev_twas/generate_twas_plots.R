library(hudson)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(RColorBrewer)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/dev_twas/rda/twas_exp_ranges.Rdata")

# Filter N/A z scores
twas_z <- twas_exp_fin %>% filter(!is.na(TWAS.Z))

twas_z_sacc <- twas_z[twas_z$region == "sacc",]

twas_z_amyg <- twas_z[twas_z$region == "amygdala",]


don <- twas_z %>%
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len = max(end)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(twas_z, ., by = c("CHR" = "CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, twas_mean_dist) %>%
    mutate(BPcum = twas_mean_dist + tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

pdf()
ggplot(don, aes(x=BPcum, y=TWAS.Z)) +

    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("#861657", "#D56AA0"), 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

    # Custom the theme:
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
dev.off()
