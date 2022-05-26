#!/bin/env Rscript

### show target insert size distribution
## load pacakges
require(ggplot2)
require(RColorBrewer)
require(plyr)

## read data
argv <- commandArgs(TRUE)
argc <- length(argv)

infile <- argv[1]
outfile <- argv[2]

insert_info <- read.delim(infile, header = TRUE, row.names = 1)

## get group mean
insert_mean <- ddply(insert_info, "insert_detect_type", summarise, grp.mean = mean(log10(insert_len)))

pdf(outfile, 8, 6)

ggplot(insert_info, aes(x = log10(insert_len), fill = insert_detect_type)) +
theme_bw() +
geom_density(alpha = 0.4) +
# ad mean lines
geom_vline(data = insert_mean, aes(xintercept = grp.mean, color = insert_detect_type), linetype = "dashed") +
scale_x_log10() +
xlab(expression(log[10]("insert size"))) +
ylab("Probability density") +
scale_fill_brewer(palette = "Set2", name = "Detection type") +
scale_color_brewer(palette = "Set2", name = "Detection type")

dev.off()
