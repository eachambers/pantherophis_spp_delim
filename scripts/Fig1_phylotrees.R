library(ggtree)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ape)
library(RColorBrewer)
library(scales)

theme_set(theme_cowplot())

## The following code generates Fig. 1, which illustrates phylo trees for both SNP and cytb datasets.
## It also generates the circular tree used in the graphical abstract/striking image for the paper.

##    FILES REQUIRED:
##          RAxML_bipartitions.Pguttatus_cytb_original.out.nexus
##          RAxML_bipartitions.result_original.nexus
##          NGSadmix/Admix_with_guttatus/k4_alldata.txt
##          NGSadmix/Admix_with_guttatus/k4_alldata_cytbannotated.txt
##          replace_labels_ngs.txt
##          replace_labels_cytb.txt
##          Code written by E. Anne Chambers

# Import trees ------------------------------------------------------------

cytb <- read.nexus("RAxML_bipartitions.Pguttatus_cytb_original.out.nexus")
ngs <- read.nexus("RAxML_bipartitions.result_original.nexus")

## Check node numbering
# ggtree(ngs) +
#   geom_text2(aes(subset = !isTip, label=node)) +
#   geom_tiplab(size=2)

## Check bootstraps
# ggtree(cytb) +
#   geom_text2(aes(subset = !isTip, label=label)) +
#   geom_tiplab(size=2)

## Build a diverging color palette for BS values, blue is high and red is low
x <- seq(0, 1, length.out = 100)
cc <- scales::seq_gradient_pal("#c7544a", "#4365d0", "Lab")(x)


# Branch (edge) colouration -----------------------------------------------

# The first nodes to be labeled are the terminal tips, which means nodes 1-135 are tips.
# You can verify this by doing the following:
# ngs_edges <- c(ngs$edge[,1], ngs$edge[,2]) # pull out nodes connecting edges
# as.data.frame(table(ngs_edges)) # count frequency of nodes; tips will only occur once

# Import NGSadmix data for K=4
ngsadmix <- read_tsv("NGSadmix/Admix_with_guttatus/k4_alldata.txt")

ngsadmix <-
  ngsadmix %>% 
  select(-ordered_number, -mtclade, -longitude) # remove irrelevant cols

# Some of the sample IDs aren't consistent between NGSadmix data and the mt tree
# so we need to import file with corrected names (FHSM and OMNH rather than F and O)
cytb_ngsadmix <- read_tsv("NGSadmix/Admix_with_guttatus/k4_alldata_cytbannotated.txt")

cytb_ngsadmix <-
  cytb_ngsadmix %>% 
  select(-ordered_number, -mtclade, -longitude)

#======== NGS TREE ===========
# Tip nodes are numbered based on the order samples occur in tree$tip.label.
# Get the node numbering for samples:
ngs_order <-
  as.data.frame(ngs$tip.label) %>% 
  rename(sampleID = "ngs$tip.label") %>% 
  mutate(node=c(1:135)) %>% 
  left_join(ngsadmix) %>% 
  gather(cluster, proportion, na.rm=FALSE, -sampleID, -node) # two samples are outgroups so have NA

# Samples with >=0.75 will be assigned the color for that cluster; others will be grey
ngs_above_threshold <-
ngs_order %>% 
  group_by(sampleID) %>% 
  filter(proportion >= 0.75) %>% 
  mutate(color = case_when(
    cluster == "k1" ~ "#be9739", # gutt
    cluster == "k2" ~ "#c7544a", # SLOW
    cluster == "k3" ~ "#75a3dd", # EMSO
    cluster == "k4" ~ "#5d8252" # EMNO
  ))

# New, brighter colors:
ngs_above_threshold <-
  ngs_order %>% 
  group_by(sampleID) %>% 
  filter(proportion >= 0.75) %>% 
  mutate(color = case_when(
    cluster == "k1" ~ "#ee0004", # gutt
    cluster == "k2" ~ "#0023f8", # SLOW
    cluster == "k3" ~ "#2d8932", # EMSO
    cluster == "k4" ~ "#efa827" # EMNO
  ))

# Figure out which samples did not have >= 0.75 proportion assignment
# For NGS data, should be 33 taxa (2 of which - KW397 and FS2838 - are outgroups)
ngs_below <- 
  as.data.frame(table(c(unique(ngs_above_threshold$sampleID), unique(ngs_order$sampleID)))) %>% 
  filter(Freq==1) %>% 
  rename(sampleID = Var1) %>% 
  mutate(color = "#969696") %>% 
  left_join(ngs_order) %>% # need to match with node number
  filter(cluster=="k1") # so there's only one value per sample

# Also add remaining (internal nodes) so they show up colored on the tree
# max(ngs$edge[,1]) # max is 269
ngs_int <- 
  as.data.frame(c(136:269)) %>% 
  rename(node = "c(136:269)") %>% 
  mutate(color = "black")

# Now, combine into a single data frame
ngs_colors <-
  full_join(ngs_above_threshold, ngs_below) %>% 
  select(-Freq) %>% 
  full_join(., ngs_int)

#======== CYTB TREE ===========
## The cytb tree has 138 tips (2 outgroups [same as NGS data] and 4 additional samples)
## These 4 additional samples are TJH1343, TJL2423, TJH3719, and TL1446.
## 2 outgroup samples (as with NGS data) are FS2838 and KW397.

# Get the node numbering for samples:
cytb_order <-
  as.data.frame(cytb$tip.label) %>% 
  rename(sampleID = "cytb$tip.label") %>% 
  mutate(node=c(1:138)) %>% 
  left_join(cytb_ngsadmix) %>% 
  gather(cluster, proportion, na.rm=FALSE, -sampleID, -node) # two samples are outgroups so have NA

# Samples with >=0.75 will be assigned the color for that cluster; others will be grey
cytb_above_threshold <-
  cytb_order %>% 
  group_by(sampleID) %>% 
  filter(proportion >= 0.75) %>% 
  mutate(color = case_when(
    cluster == "k1" ~ "#be9739", # gutt
    cluster == "k2" ~ "#c7544a", # SLOW
    cluster == "k3" ~ "#75a3dd", # EMSO
    cluster == "k4" ~ "#5d8252" # EMNO
  ))

# Figure out which samples did not have >= 0.6 proportion assignment
# For NGS data, should be 37 taxa (2 of which are outgroups)
cytb_below <- 
  as.data.frame(table(c(unique(cytb_above_threshold$sampleID), unique(cytb_order$sampleID)))) %>% 
  filter(Freq==1) %>% 
  rename(sampleID = Var1) %>% 
  mutate(color = "grey") %>% 
  left_join(cytb_order) %>% # need to match with node number
  filter(cluster=="k1") # so there's only one value per sample

# Also add remaining (internal nodes) so they show up colored on the tree
# max(ngs$edge[,1]) # max is 269
cytb_int <- 
  as.data.frame(c(139:275)) %>% 
  rename(node = "c(139:275)") %>% 
  mutate(color = "black")

# Now, combine into a single data frame
cytb_colors <-
  full_join(cytb_above_threshold, cytb_below) %>% 
  select(-Freq) %>% 
  full_join(., cytb_int)

# Build NGS tree (right-handed, default) ----------------------------------------------------------

# First, color branches according to ngs_colors
ngs_cols <-
  ggtree(ngs) %<+%
  ngs_colors + aes(color=I(color)) +
  geom_tiplab(size=2)

ngs_cols <-
  ngs_cols %>% 
  flip(158, 171)

# Now, build tree with aesthetics we want
ngs_tree <-
  ggtree(ngs) +
  coord_cartesian(xlim = c(0,.25), clip="off") +
  # geom_tiplab(size=2) +
  # geom_hilight(node=138, fill='#be9739') + # GUTT
  # geom_hilight(node=148, fill='#c7544a') + # SLOW
  # geom_hilight(node=208, fill='#5785bf') + # EMSO
  # geom_hilight(node=269, fill='#b2b2b2') + # OBSO
  # geom_hilight(node=233, fill='#5d8252') + # EMWE; same colour as EMNO
  # geom_hilight(node=242, fill='#5d8252') + # EMNO
  # geom_text2(aes(subset = !isTip, label=label)) +
  geom_nodepoint(aes(color=as.numeric(label), size=as.numeric(label))) +
  scale_color_gradientn(colors = cc) +
  geom_treescale() +
  scale_size(range = c(0.1, 2))
  
ngs_tree <-
  ngs_tree %>% 
  flip(158, 171)

## CHANGE SOME OF THE LABELS
# Can take a look at existing labels on ngs tree & their ordering
ggtree:::get.tree(ngs_tree)$tip.label

# Import replacement labels with corrected names
ngs_labels <- read_tsv("replace_labels_ngs.txt")

# Build tree with corrected labels (ensure original tree doesn't have labels at all)
ngs_tree <-
  ngs_tree %<+% ngs_labels + geom_tiplab(aes(label=new_label), size=2)

# Build left-handed cytb tree ----------------------------------------------

# To get the cytb tree in the lefthand direction, we need to adjust the axes values
# First make the default cytb tree
cytb_tree <- ggtree(cytb)

# Extract data from default cytb tree & assign them to a new dataframe
rev_data <- cytb_tree$data
rev_data$x <- max(rev_data$x) - rev_data$x # reverse the values

# ggtree(rev_data) +
#   geom_text2(aes(subset = !isTip, label=label)) +
#   geom_tiplab(size=2, hjust = 1.3)

# First, color branches according to cytb_colors
cytb_cols <-
  ggtree(rev_data) %<+%
  cytb_colors + aes(color=I(color)) +
  geom_tiplab(size=2, hjust = 1.3)

cytb_cols %>% 
flip(253, 194) %>% 
  flip(211, 195)

#===============
new_cytb <-
ggtree(rev_data) +
  coord_cartesian(xlim = c(0,1), clip="off") +
  # geom_tiplab(size=1.75, hjust=1.3) +
  # geom_hilight(node=185, fill='#be9739') + # GUTT
  # geom_hilight(node=142, fill='#c7544a') + # SLOW
  # geom_hilight(node=253, fill='#5785bf') + # EMSO
  # geom_hilight(node=275, fill='#b2b2b2') + # OBSO
  # geom_hilight(node=195, fill='#5d8252') + # EMWE
  # geom_hilight(node=211, fill='#5d8252') + # EMNO
  geom_nodepoint(aes(color=as.numeric(label), size=as.numeric(label))) +
  scale_colour_gradientn(colors=cc) +
  geom_treescale() +
  scale_size(range = c(0.1, 2))

new_cytb <-
  new_cytb %>% 
  flip(253, 194) %>% 
  flip(211, 195)

## CHANGE SOME OF THE LABELS
# Can take a look at existing labels on ngs tree & their ordering
as.data.frame(ggtree:::get.tree(new_cytb)$tip.label)

# Import replacement labels with corrected names
cytb_labels <- read_tsv("replace_labels_cytb.txt")

# Build tree with corrected labels (ensure original tree doesn't have labels at all)
new_cytb <-
  new_cytb %<+% cytb_labels + geom_tiplab(aes(label=new_label), size=1.75, hjust=1.3)

# Export trees ------------------------------------------------------------

# ngs_cols
# ggsave("Pantherophis_NGS_tree_coloredbranches.pdf", width=7.8, height = 8.3)
ngs_tree
ggsave("Pantherophis_NGS_tree_NOTIPS.pdf", width=7.8, height = 8.3)

# cytb_cols
# ggsave("Pantherophis_CYTB_tree_coloredbranches.pdf", width=7.8, height = 8.3)
new_cytb
ggsave("Pantherophis_CYTB_tree.pdf", width=7.8, height = 8.3)


# STRIKING IMAGE ----------------------------------------------------------
circ <-
ggtree(ngs, layout='circular', size=1) %<+%
  ngs_colors + aes(color=I(color)) 

br_colored <-
ggtree(ngs, size=1) %<+%
  ngs_colors + aes(color=I(color)) 

circ
ggsave("Striking_tree_circular.pdf", width=7.8, height = 8.3)
br_colored
ggsave("Striking_tree_brcolored.pdf", width=7.8, height = 8.3)

