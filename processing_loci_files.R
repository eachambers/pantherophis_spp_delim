setwd("~/Box Sync/Lampropeltis project/PANTHEROPHIS_project/IQtree_PHRAPL")

library(tidyverse)
library(cowplot)
library(ggplot2)
library(phylotools)

## The following code takes in files for each locus (generated using the Split_loci_file.ipynb
## script) and makes them into phylip files. It then randomly selects loci for downstream tree
## reconstruction for PHRAPL.
##  Written by E. Anne Chambers.

##    (a) Get .loci file output from iPyrad.
##    (b) Run Split_loci_file.ipynb and save separate loci files in pantherophis_loci_output.
##  > (c) Run processing_loci_files.R to get separate phylip files for each locus and randomly subsample files.
##    (d) Run count_pop_loci.R to count inds per population.
##    (e) Pull out relevant models from migrationArray and save as new migrationArray.


##    FILES REQUIRED:
##          pantherophis_loci_output/pantherophis_locus*.phy


# CONVERT TO PHYLIP -------------------------------------------------------

## Function to read in and convert file to phylip
convert_file <- function (filename) { # can do browser() to run through each line of function
  # read in files and make temp object t
  t <- read.table(filename, header=F, sep="")
  # want to change output files to be diff from original files
  outfile_name <- str_replace(string = filename, pattern = ".phy", replacement = ".phylip")
  # convert to phylip and save as new output file names
  dat2phylip(t, outfile=outfile_name)
}

## Test function on a single file
# convert_file("pantherophis_loci_output_test/pantherophis_locus28.phy")

## List of all files to do function on
file_list <- list.files(path="pantherophis_loci_output/", pattern="*.phy", full.names=TRUE, recursive=FALSE) # 7021 elements (# loci)

## Run function on list of files
map(file_list, convert_file) # new version of lapply; no return statement so output is null


# TO RANDOMLY SUBSAMPLE LOCI -----------------------------------------------------

# From Bruno Silva (https://www.r-bloggers.com/2019/04/random-sampling-of-files/)
# will move randomly selected files (along with a .csv with file list) to selected/
# directory within chosen path

random_files <- function(path, percent_number, pattern = "phylip"){
  ####################################################################
  # path = path to folder with files to select                                 
  #                                                                            
  # percent_number = percentage or number of recordings to select. If value is 
  #   between 0 and 1 percentage of files is assumed, if value greater than 1, 
  #   number of files is assumed                                               
  #                                                                            
  # pattern = file extension to select      
  ####################################################################
  
  # Get file list with full path and file names
  files <- list.files(path, full.names = TRUE, pattern = pattern)
  file_names <- list.files(path, pattern = pattern)
  
  # Select the desired % or number of file by simple random sampling 
  randomize <- sample(seq(files))
  files2analyse <- files[randomize]
  names2analyse <- file_names[randomize]
  if(percent_number <= 1){
    size <- floor(percent_number * length(files))
  }else{
    size <- percent_number
  }
  files2analyse <- files2analyse[(1:size)]
  names2analyse <- names2analyse[(1:size)]
  # Create folder to output
  results_folder <- paste0(path, '/selected')
  dir.create(results_folder, recursive=TRUE)
  
  # Write csv with file names
  write.table(names2analyse, file = paste0(results_folder, "/selected_files.csv"),
              col.names = "Files", row.names = FALSE)
  
  # Move files
  for(i in seq(files2analyse)){
    file.rename(from = files2analyse[i], to = paste0(results_folder, "/", names2analyse[i]) )
  }
}

# Run function for 200 loci files:
random_files(path ="pantherophis_loci_output/", percent_number = 200, pattern = "phylip")
