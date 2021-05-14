setwd("~/Box Sync/Lampropeltis project/PANTHEROPHIS_project")

library(scatterpie) # needed for making pies
library(cowplot)
library(tidyverse)

theme_set(theme_cowplot())

## The following code generates Fig. S2, which illustrates piecharts indicating ancestry 
##  from NGSadmix analysis
##  Written by E. Anne Chambers

##    FILES REQUIRED:
##          Map_q_input.csv

# Import data -------------------------------------------------------------

# Upload data and define variables
dat <- read.csv("Map_q_input.csv")

us <- map_data("state")

# Build main (large) plot ---------------------------------------------------------

## q1 = guttatus (ochre) "#be9739"
## q2 = slowinskii (red) "#c7544a"
## q3 = emoryi south (blue) "#75a3dd"
## q4 = emoryi north (green) "#5d8252"

# New, brighter colors
## q1 = guttatus (red) "#ee0004"
## q2 = slowinskii (blue) "#0023f8"
## q3 = emoryi south (green) "#2d8932"
## q4 = emoryi north (yellow) "#efa827"

# Pies need to be manually jittered after the fact; for these points,
# we'll need the points plotted (in black) as well as the pies
lrg_overlaps <- c("TJL2431",
                  "TJL1872",
                  "TJL3121",
                  "TJL2893",
                  "TJL2894",
                  "TJL2217",
                  "TJL3029",
                  "TJL3028",
                  "TJL3083",
                  "O43445",
                  "O43444",
                  "O43943",
                  "O43819",
                  "O42942",
                  "O42939",
                  "O41875",
                  "O42946",
                  "O42945",
                  "F17580",
                  "F17579",
                  "CSA179",
                  "CSA167",
                  "MTH440",
                  "UTA58667",
                  "H15909",
                  "H19690",
                  "H16057",
                  "H20638",
                  "H8537")

lrg_jitter <-
  dat %>% 
  filter(ind %in% lrg_overlaps) # should have 29 obs

# Build the background map
plot <-
ggplot() + 
  geom_polygon(data=us, aes(x=long, y=lat, group=group),
                color="grey", fill="#e6e6e6", size=1) +
  coord_equal(xlim = c(-108.5,-81),ylim = c(25,40)) +
  # geom_point(data = dat, aes(x=lon, y=lat)) +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank())

# Build the main map plot with certain points in black that'll be jittered later
main_plot <-
  plot + 
  geom_scatterpie(data=dat, aes(x=lon, y=lat, group=ind),
                    cols=2:5, pie_scale = 0.7, color=NA, alpha=0.9) + # alpha=0.8 for duller colors 
  theme(legend.position="none") +
  # scale_fill_manual(values = c("q1"="#be9739", "q2"="#c7544a", "q3"="#75a3dd", "q4"="#5d8252")) +
  scale_fill_manual(values = c("q1"="#ee0004", "q2"="#0023f8", "q3"="#2d8932", "q4"="#efa827")) + # brighter colors
  geom_point(data = lrg_jitter, aes(x=lon, y=lat), size=3) +
  geom_rect(aes(xmin=-101, xmax=-96, ymin=26, ymax=31), color="grey", fill=NA, size=1)

# Build zoomed in (S Texas) panel ----------------------------------------------------

# These are data within a certain geographic range
jitter <- 
dat %>% 
  filter(lon >= -98.1 & lon <= -96) %>% 
  filter(lat >= 29.6 & lat <= 31)

# Add on a couple other points that have overlapping pies
overlaps <- c("FS3061",
              "CSA272",
              "CSA563",
              "TJH3395",
              "TJH3394",
              "CSA485",
              "CSA472",
              "TM166",
              "TM165",
              "TJL3082",
              "DRD5478",
              "DRD5477",
              "TM154",
              "DRD068",
              "TJL3005")

more_jitter <-
dat %>% 
  filter(ind %in% overlaps)

# Combine the two together
jitter <-
  full_join(jitter, more_jitter) # should be 50 data points

## Build plots themselves
zoom <-
  ggplot() + 
  geom_polygon(data=us, aes(x=long, y=lat, group=group),
               color="grey", fill="#e6e6e6" ) +
  coord_equal(xlim = c(-101,-96),ylim = c(26,31)) +
  theme(
    # plot.background = element_rect(colour="grey", size=1),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank())

zoomed <-
zoom + 
  geom_scatterpie(data=dat, aes(x=lon, y=lat, group=ind),
                  cols=2:5, color=NA, pie_scale = 0.3, alpha=0.8) +
  theme(legend.position="none") +
  # scale_fill_manual(values = c("q1"="#be9739", "q2"="#c7544a", "q3"="#75a3dd", "q4"="#5d8252")) +
  scale_fill_manual(values = c("q1"="#ee0004", "q2"="#0023f8", "q3"="#2d8932", "q4"="#efa827")) + # brighter colors
  geom_point(data = jitter, aes(x=lon, y=lat), size=3)
  # geom_text(data=dat, aes(x=lon, y=lat, label=ind), size=2)

# Export plots ------------------------------------------------------------

main_plot
ggsave("FigS2_main_figure.pdf", width=16.625, height = 9.24)

zoomed
ggsave("FigS2_inset_panel.pdf", width=7.030, height = 6.798)
