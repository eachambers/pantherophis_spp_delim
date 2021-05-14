setwd("~/NGSadmix_files/")

library(pophelper)
library(gridExtra)
library(cowplot)
library(tidyverse)

theme_set(theme_cowplot())

## The following code generates Fig. S1, which illustrates NGSadmix results as Structure-style
##  plots. Code written by E. Anne Chambers

##    FILES REQUIRED:
##          KX_alldata.txt - these are the results for each K value, sorted by mitochondrial clade and then longitude
##          KX_sorted.txt - these are the results from the subsampled dataset (Fig. S1, panel B)
##          (Original order of inds and membership is in the inds2pops and inds2pops_sorted files)


# Import data -------------------------------------------------------------

# All data (Fig. S1, panel A)
k2_data <- read_tsv("k2_alldata.txt", col_names = TRUE)
k3_data <- read_tsv("k3_alldata.txt", col_names = TRUE)
k4_data <- read_tsv("k4_alldata.txt", col_names = TRUE)
k5_data <- read_tsv("k5_alldata.txt", col_names = TRUE)

# Subsampled data (Fig. S1, panel B)
k2_sub <- read_tsv("k2_sorted.txt", col_names=TRUE)
k3_sub <- read_tsv("k3_sorted.txt", col_names=TRUE)
k4_sub <- read_tsv("k4_sorted.txt", col_names=TRUE)

# Tidy data ---------------------------------------------------------------

tidy_structure_data <- function(data){
  data %>% 
    gather(key="cluster", value="proportion", -ordered_number, -sampleID, -mtclade, -longitude)
}

k2_data <- tidy_structure_data(k2_data)
k3_data <- tidy_structure_data(k3_data)
k4_data <- tidy_structure_data(k4_data)
k5_data <- tidy_structure_data(k5_data)

tidy_structure_data_sub <- function(data){
  data %>% 
    gather(key="cluster", value="proportion", -ordered_number, -sampleID, -mtclade)
}

k2_sub_new <- tidy_structure_data_sub(k2_sub)
k3_sub_new <- tidy_structure_data_sub(k3_sub)
k4_sub_new <- tidy_structure_data_sub(k4_sub)

## Assign colors to clusters
k2_cols=c("#ee0004", "#2d8932")
k3_cols=c("#efa827", "#ee0004", "#2d8932")
k4_cols=c("#ee0004", "#0023f8", "#2d8932", "#efa827")
k5_cols=c("#2d8932", "#efa827", "#0023f8", "#c7d5a9","#ee0004")


# Build the plot for panel A ----------------------------------------------------------

structure_plot <- function(data, colors) {
  data %>% 
  ggplot(aes(x=ordered_number, y=proportion, fill=cluster)) +
    geom_bar(stat="identity") + 
    panel_border() + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_manual(values=colors) +
    theme_cowplot() %+replace% theme(axis.line=element_line(colour="black"),
                                     axis.text.x=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.ticks=element_blank(),
                                     axis.title.x=element_blank(),
                                     axis.title.y=element_blank(),
                                     legend.position="none",
                                     # panel.border = element_rect(fill=NA, colour="black", linetype="solid",size=1.5),
                                     strip.text.y = element_text(size=30, face="bold"),
                                     strip.background = element_rect(colour="white", fill="white"),
                                     panel.spacing=unit(-0.1, "lines")) +
    geom_vline(xintercept = 23.5, size=1) +
    geom_vline(xintercept = 81.5, size=1) +
    geom_vline(xintercept = 124.5, size=1)
}


p_k2 <- structure_plot(k2_data, k2_cols)
p_k3 <- structure_plot(k3_data, k3_cols)
p_k4 <- structure_plot(k4_data, k4_cols)
p_k5 <- structure_plot(k5_data, k5_cols)

full_plot <- plot_grid(p_k2, p_k3, p_k4, p_k5, nrow = 4)
full_plot
ggsave("FigS1_panelA_NGSadmix.pdf", width=8.6, height = 7.856)


# Build the plot for panel B ----------------------------------------------

structure_plot_sub <- function(data, colors) {
  data %>% 
    ggplot(aes(x=ordered_number, y=proportion, fill=cluster)) +
    geom_bar(stat="identity") + 
    panel_border() + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_manual(values=colors) +
    theme_cowplot() %+replace% theme(axis.line=element_line(colour="black"),
                                     axis.text.x=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.ticks=element_blank(),
                                     axis.title.x=element_blank(),
                                     axis.title.y=element_blank(),
                                     legend.position="none",
                                     # panel.border = element_rect(fill=NA, colour="black", linetype="solid",size=1.5),
                                     strip.text.y = element_text(size=30, face="bold"),
                                     strip.background = element_rect(colour="white", fill="white"),
                                     panel.spacing=unit(-0.1, "lines")) +
    geom_vline(xintercept = 9.5, size=1) +
    geom_vline(xintercept = 18.5, size=1) +
    geom_vline(xintercept = 27.5, size=1)
}


p_k2_sub <- structure_plot_sub(k2_sub_new, k2_cols)
p_k3_sub <- structure_plot_sub(k3_sub_new, k3_cols)
p_k4_sub <- structure_plot_sub(k4_sub_new, k4_cols)

full_plot_sub <- plot_grid(p_k2_sub, p_k3_sub, p_k4_sub, nrow = 3)
full_plot_sub
ggsave("FigS1_panelB_NGSadmix.pdf", width=8.6, height = 7.856)

