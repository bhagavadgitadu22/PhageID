#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: plot_long_viral_reads.R <geNomad summary.tsv> <plot.png>")
}

suppressPackageStartupMessages(library(ggplot2))
histogram <- read.delim(args[1], check.names = FALSE)
if (!"length" %in% names(histogram)) {
  stop("geNomad summary must contain a length column")
}
histogram$length <- suppressWarnings(as.numeric(histogram$length))
if (anyNA(histogram$length)) {
  stop("geNomad summary must contain a numeric length column")
}

if (nrow(histogram) == 0) {
  plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No viral reads >=20 kbp") +
    theme_void()
} else {
  plot <- ggplot(histogram, aes(x = length)) +
    geom_histogram(binwidth = 100) +
    labs(x = "Length (bp)", y = "Count") +
    theme_classic()
}

dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
ggsave(args[2], plot = plot, width = 8, height = 4.5, dpi = 300, bg = "white")
