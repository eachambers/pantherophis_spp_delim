library(devtools)
library(ape)
library(tidyverse)
library(phangorn)
library(phrapl)

# devtools::install_github("bomeara/phrapl")

## The following code does two processing steps, and then runs PHRAPL. It roughly follows the format of Tutorial 1
## within the PHRAPL documentation (https://phrapl.github.io/Content/Tutorials/ModelSelection/6aa.PHRAPL.tutorial1.html).
##  Written by E. Anne Chambers.

## Briefly, it:
##    (a) Removes outgroups and branch lengths from trees and save.
##    (b) Calculates the number of individuals per population per tree (locus).
##    (c) Removes irrelevant trees.
##    (d) Gets subsampling weights from trees, builds observedTrees and saves PHRAPL input file.


##    FILES REQUIRED:
##          raw output best trees from RAxML in a directory (not provided in supp materials)
##          master tree file (done in bash after step a complete): pantherophis_master_trees.tre
##          population assignment file: pantherophis_N4_popass.txt
##          migration array: MigrationArray_4Pop_migonlyadjacent_2944models_K1.rda


# a. REMOVE BLS AND OGS FROM TREES ---------------------------------------------------

## Function to remove branch lengths and outgroups from tree files
remove_bls_ogs <- function (filename) {
  # browser()
  # read in files and make temp object t
  t <- read.tree(filename)
  # remove outgroup samples (if present)
  t <- drop.tip(t, c("KW397", "FS2838"))
  # remove branch lengths
  t$edge.length <- NULL
  t$root.edge <- NULL
  # want to change output files to be diff from original files
  outfile_name <- str_replace(string = filename, pattern = ".phylip.raxml.bestTree", replacement = "_noogs_nobls.tre")
  # convert to phylip and save as new output file names
  write.tree(t, file=outfile_name)
}

## Test function on a single file
# remove_bls_ogs("best_trees_selected/pantherophis_locus424.phylip.raxml.bestTree")

## List of all files to do function on and run function on list
file_list <- list.files(path="all_trees", pattern="*.bestTree", full.names=TRUE, recursive=FALSE) # 6249 elements (# trees)
map(file_list, remove_bls_ogs)

## Now, use bash in terminal to combine trees into a single file (faster if done in terminal:
# system('cat best_trees_selected/*_noogs_nobls.tre >> pantherophis_master_trees.tre')
# system('cat all_trees/*_noogs_nobls.tre >> pantherophis_master_trees.tre')

# b. CALCULATE NO. INDS PER POPULATION PER TREE -------------------------------------

currentAssign <- read.table("pantherophis_N4_popass.txt", sep="\t", header=TRUE) # 133 obs
currentTrees <- read.tree("pantherophis_master_trees.tre") # 6445 elements

# We can see tip labels by doing this:
# currentTrees[[1]][3] # or
# currentTrees[[1]]$tip.label

# Need to iterate through all trees within currentTrees and match against the assignment file
## Function to count tips and match to populations from each tree
count_phylo_pop <- function (phyloobject, treeid, assignment=currentAssign) {
  # browser()
  if(is.null(phyloobject$tip.label)) browser()
  tibble(sampleID = phyloobject$tip.label) %>% 
    inner_join(assignment, by = 'sampleID') %>% 
    count(poplabel) %>% 
    mutate(tree = treeid)
}

# Check that this works:
# count_phylo_pop(currentTrees[[1]], treeid = 1, assignment = currentAssign)

# use map2 because there are two inputs that are changing
# map2 will automatically go through lists x (currentTrees) and y (seq_along) and
# run whatever function along; can also ... after function specified if you want
# to add assignment = currentAssign so it applies it every time it maps, but 
# because it's specified in the function itself it's unnecessary
matched_trees_taxa <- map2(currentTrees, 
                            seq_along(currentTrees), 
                            count_phylo_pop) %>% 
  bind_rows() # 25777 obs (would be 6445*4=25780 if all loci had 4 populations)

# matched_trees_taxa is a dataframe indicating each tree, the population, and the number
# of individuals within that population.


# c. REMOVE IRRELEVANT TREES -----------------------------------------------

# Trees without individuals from a given population must be removed.
# Only necessary if you haven't run count_pop_loci.R on phylip files before reconstructing trees.
# In our case, only ~15 trees need to be removed from that original 196 locus dataset

no_NAs_taxa <- # now 25768 obs (6442 trees, 3 fewer than original)
  matched_trees_taxa %>% 
  pivot_wider(names_from = poplabel, values_from = n) %>% 
  drop_na() %>% # remove rows (=trees) with NAs (i.e., no inds in a pop'n)
  pivot_longer(names_to="poplabel", values_to="n", cols=2:5)

# Now, let's find which trees have fewer than a certain number of specimens for any population.
min_no_taxa <- 5
threshold_taxa <- # 12 trees
  matched_trees_taxa %>% 
  filter(n<min_no_taxa)

no_NAs_taxa <-
  no_NAs_taxa %>% 
    filter(!(no_NAs_taxa$tree %in% threshold_taxa$tree)) # 25720 obs; 6430 trees

# Remove trees
new_currentTrees <-
  currentTrees[unique(no_NAs_taxa$tree)] # 6430 trees; 15 fewer than we began with

# Subsamples trees to only 100 loci
random100 <- sample(c(1:6430), 100)

new_currentTrees[[]]

# d. SUBSAMPLING ----------------------------------------------------------

##  Much of the following PHRAPL code is from Ariadna Morales from 
##  Morales et al. 2017 Syst. Biol. 66(3):440-452. doi:10.1093/sysbio/syw100

## Arguments ##
nloci <- length(new_currentTrees)     ## Define number of loci.
popAssignments <- list(c(2,2,2,2))    ## Subsamples per pop; sims. showed 2-3 inds per pop is enough.
subsamplesPerGene <- min_no_taxa      ## Subsamples per gene; the number of replicate subsamples to take per locus

## Load input data ##
# Assignment file
currentAssign <- read.table("pantherophis_N4_popass.txt", sep="\t", header=TRUE) # 133 obs

## Do subsampling ##
observedTrees <- PrepSubsampling(assignmentsGlobal=currentAssign,
                               observedTrees=new_currentTrees,
                               popAssignments=popAssignments,
                               subsamplesPerGene=subsamplesPerGene,
                               outgroup=FALSE,
                               outgroupPrune=FALSE) # 905 trees

# observedTrees is a list of subsampled trees that match format of popAssignments
# (one set per popVector)
# to access first subsample iteration of second locus:
# observedTrees[[1]][[11]] # 6 tips and 5 internal nodes

## Get subsample weights ##
subsampleWeights.df <- GetPermutationWeightsAcrossSubsamples(popAssignments=popAssignments,
                                                           observedTrees=observedTrees) # 32150 trees

## Save RDA file ##
save(list=c("observedTrees","subsampleWeights.df"),file="pantherophis_phraplInput_4pops_alltrees.rda")