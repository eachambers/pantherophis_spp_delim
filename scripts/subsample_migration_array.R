library(devtools)
library(ape)
library(tidyverse)
library(phangorn)
library(phrapl)

##  Written by E. Anne Chambers.
## The following code takes in a migration array and pulls out models within it relevant to our analysis.
## In our case, we only want to consider models that have migration between geographically adjacent populations.
## This means we're going to allow migration between the following:
##    A(EMSO)<->B(EMNO)
##    A(EMSO)<->C(SLOW)
##    B(EMNO)<->C(SLOW)
##    C(SLOW)<->D(GUTT)

##    (a) Get .loci file output from iPyrad.
##    (b) Run Split_loci_file.ipynb and save separate loci files in pantherophis_loci_output.
##    (c) Run processing_loci_files.R to get separate phylip files for each locus and randomly subsample files.
##    (d) Run count_pop_loci.R to count inds per population.
##  > (e) Pull out relevant models from migrationArray and save as new migrationArray.
##    (f) Remove outgroups and branch lengths from trees and save.
##    (g) Calculate the number of individuals per population per tree (locus).


##    FILES REQUIRED:
##          migration array: MigrationArray_4pops_maxMigrationK1; this is pre-loaded as part of PHRAPL package


# LOAD MIGRATION ARRAY ----------------------------------------------------

# See which migration arrays are available:
# unlist(strsplit(grep("MigrationArray",list.files(paste(path.package("phrapl"),  
#                                                        "/data/",sep="")),value=TRUE),".rda"))

# Load one of these (it'll be saved as migrationArray):
data("MigrationArray_4pops_maxMigrationK1") # equal pop'n size, resolved trees, symm migration, 2 migration rates
# length(migrationArray) # 18,432 models

# Population assignment is as follows: A=EMSO; B=EMNO; C=SLOW; D=GUTT
# Allow migration between spatially adjacent pop'ns:
# A/B (EMSO/EMNO);  A/C (EMSO/SLOW);  B/C (EMNO/SLOW);  C/D (SLOW/GUTT)

# Cells that NEED to have 0 (we cannot possibly have migration)
na_cells <- c(4, 8, 13, 14, 20, 24, 29, 30, 36, 40, 45, 46)

## Function to determine if na_cells has NAs or 0 (only models we're interested in)
check_mig_nas <- function (array, na_cells) { # can do browser() to run through each line of function
  # browser()
  all(is.na(array$migrationArray[na_cells]) | array$migrationArray[na_cells]==0)
}

# Check that this works:
# check_mig_nas(migrationArray[[2]], model=2)

# Run on all models within migrationArray
only_nas <-
  tibble(ma = migrationArray) %>% 
  mutate(index = seq_along(ma),
         is_relevant = map(ma, check_mig_nas, na_cells=na_cells) %>% unlist()) %>% 
  filter(is_relevant==TRUE)

# Get a count of each
# only_nas %>% count(is_relevant)


# EXTRACT RELEVANT MODELS AND SAVE -------------------------------------------------

MigrationArray <- only_nas$ma
MigrationArrayMap <- GenerateMigrationArrayMap(MigrationArray) 
save(MigrationArray, MigrationArrayMap, file="MigrationArray_4Pop_migonlyadjacent_2944models_K1.rda")

# Check that this loads properly
# load("MigrationArray_4Pop_migonlyadjacent_2944models_K1.rda")


# SUBSAMPLE FURTHER -------------------------------------------------------

# Let's now just take 20 random models from each (of 18) topologies (MigrationArrayMap has this info in the
# collapseMatrix.number column)

load("MigrationArray_4Pop_migonlyadjacent_2944models_K1.rda")

# randomly sample 20 models from each topology (collapseMatrix.number 1-18)
# arranged in dataframe such that one model per topology will run at a time (rather than all models of one topology)
random_subsamples <- 
  MigrationArrayMap %>% 
  group_by(collapseMatrix.number) %>% 
  sample_n(size=20) %>% 
  mutate(index = seq_along(model)) %>% 
  arrange(index) %>% 
  select(model, collapseMatrix.number, index) # randomly arrange models in job script

random_models <-
  random_subsamples$model


# SAMPLES FOR PERSONAL COMP -----------------------------------------------

load("MigrationArray_4Pop_migonlyadjacent_2944models_K1.rda")

# randomly sample 5 models from topologies 9-18 (collapseMatrix.number 9-18)
fewer_subsamples <- 
  MigrationArrayMap %>% 
  group_by(collapseMatrix.number) %>% 
  filter(collapseMatrix.number == c(9:18)) %>% 
  sample_n(size=5) %>% # 50 observations
  filter(model != 587, model != 597, model != 712) %>% # should be 47 obs
  filter(collapseMatrix.number != 10) # 42 obs

fewer_models <-
  fewer_subsamples$model

