library(devtools)
library(ape)
library(tidyverse)
library(phangorn)
library(purrr)


## The following code takes in phylip files for each locus and counts how many individuals per
## population there are in each. This is required because PHRAPL cannot run if some loci don't 
## have all/adequate numbers of individuals per population in each locus.
##  Written by E. Anne Chambers.

##    (a) Get .loci file output from iPyrad.
##    (b) Run Split_loci_file.ipynb and save separate loci files in pantherophis_loci_output.
##    (c) Run processing_loci_files.R to get separate phylip files for each locus and randomly subsample files.
##  > (d) Run count_pop_loci.R to count inds per population.
##    (e) Pull out relevant models from migrationArray and save as new migrationArray.
##    (f) Remove outgroups and branch lengths from trees and save.
##    (g) Calculate the number of individuals per population per tree (locus).


##    FILES REQUIRED:
##          population assignment file: pantherophis_N4_popass.txt
##          path to all loci (phylip) files you want counted: e.g., pantherophis_loci_output/other_files (not provided)
##          path to loci you've already counted: e.g., pantherophis_loci_output/selected/ (not provided)

# CALCULATING NO. INDS PER POPULATION PER TREE -------------------------------------

currentAssign <- read.table("pantherophis_N4_popass.txt", sep="\t", header=TRUE) # 133 obs

# Need to iterate through all trees within currentTrees and match against the assignment file
## Function to count tips and match to populations from each tree
count_pop_loci <- function (filename, assignment=currentAssign) { # can do browser() to run through each line of function
  # browser()
  t <- read.table(filename, sep="" , header=F)
  t = t[-1,] # delete first row of file which contains no. taxa and no. sites
  # remove outgroup samples
  t <-
    t %>% filter(V1 != c("KW397", "FS2838"))
  # now compare t to assignment file
  t %>% 
    rename(sampleID = V1) %>% 
    inner_join(assignment, by = 'sampleID') %>% 
    count(poplabel) %>% 
    mutate(locus = filename)
}

# Check that this works:
# count_pop_loci(filename="pantherophis_loci_output/other_files/pantherophis_locus6794.phylip")
## ^ should have A=19, B=49, C=34, and D=7

## List of all files to do function on
## Be sure there aren't any directories within this one with the common string "phy"
file_list <- list.files(path="pantherophis_loci_output/other_files", 
                        pattern="*phylip", 
                        full.names=TRUE, 
                        recursive=FALSE) # 7021 elements (# loci)

## Run function on list of files
no_taxa_loci <-
  map(file_list, count_pop_loci) %>% 
  bind_rows() # error messages will be spit out; this is fine- just happens when a locus doesn't have any outgroup samples
# ^ 27885 obs; 7021*4=28084 so there are ~199 loci without all four populations


# GET STATS ---------------------------------------------------------------

# Loci without individuals from a given population
loci_with_all_pops <- # now 27288 obs (6822 loci, 199 fewer than original)
  no_taxa_loci %>% 
  pivot_wider(names_from = poplabel, values_from = n) %>% 
  drop_na() %>% # remove rows (=loci) with NAs (i.e., no inds in a pop'n)
  pivot_longer(names_to="poplabel", values_to="n", cols=2:5)

# Now, let's take a look at which loci are beneath the threshold for the min_no_taxa
min_no_taxa <- 5

threshold_loci <- # 373 loci (so total will be 6822-373=6449 loci)
  loci_with_all_pops %>% 
  filter(n<min_no_taxa)

loci_above_threshold <-
  loci_with_all_pops %>% 
  filter(!(loci_with_all_pops$locus %in% threshold_loci$locus)) # 25796 obs; 25796/4=6449 loci

# Remove directory information
loci_above_threshold$locus <-
  gsub("pantherophis_loci_output/other_files/", "", loci_above_threshold$locus)

# Export relevant column
write.table(unique(loci_above_threshold$locus), "pantherophis_loci_output/other_files/loci_above_threshold.txt", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)



# DON'T REBUILD TREES  ----------------------------------------------------

# We've already constructed trees for some loci; we don't want to waste time doing so again.
# We need to match loci for trees we've built with those in loci_above_threshold.
# List of loci we've already built trees for
file_list_done <- list.files(path="pantherophis_loci_output/selected", 
                             pattern="*phylip", 
                             full.names=TRUE, 
                             recursive=FALSE) # 200 elements (# trees)

no_taxa_loci_done <-
  map(file_list_done, count_pop_loci) %>% 
  bind_rows() # 796 obs; 200*4=800 so there are ~4 loci without all four populations

loci_with_all_pops_done <- # now 784 obs (196 loci, 4 fewer than original)
  no_taxa_loci_done %>% 
  pivot_wider(names_from = poplabel, values_from = n) %>% 
  drop_na() %>% # remove rows (=loci) with NAs (i.e., no inds in a pop'n)
  pivot_longer(names_to="poplabel", values_to="n", cols=2:5) # 784 obs

threshold_loci_done <- # 13 loci (so total will be 196-13=183 loci)
  loci_with_all_pops_done %>% 
  filter(n<min_no_taxa)

loci_above_threshold_done <-
  loci_with_all_pops_done %>% 
  filter(!(loci_with_all_pops_done$locus %in% threshold_loci_done$locus)) # 732 obs; 732/4=183 loci

# Remove directory information
loci_above_threshold_done$locus <-
  gsub("pantherophis_loci_output/selected/", "", loci_above_threshold_done$locus)


# REMOVE ROWS OF LOCI THAT HAVE TREES -------------------------------------

loci_above_threshold_new <-
  loci_above_threshold %>%
    filter(!(loci_above_threshold$locus %in% loci_above_threshold_done$locus)) # 25796-732=25064 obs.; 25064/4=6266 loci

# Export relevant column
write.table(unique(loci_above_threshold_new$locus), "pantherophis_loci_output/other_files/without_finished_trees.txt", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)
