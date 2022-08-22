# pantherophis_spp_delim

The following is code to generate figures from [Marshall et al. (2021)](https://www.sciencedirect.com/science/article/abs/pii/S1055790321001275). Raw data files are provided on Dryad [here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.cc2fqz64d).

## Data analysis
* [Split iPyrad .loci files and save separately](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/Split_loci_file.ipynb)
* [Generate .phylip file for each locus for input into RAxML-ng](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/processing_loci_files.R)
* [Count number of individuals per locus for input into PHRAPL](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/count_pop_loci.R)
* [Subsample pre-loaded PHRAPL MigrationArray](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/subsample_migration_array.R)
* [Additional data processing and generate PHRAPL input file](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/PHRAPL_workthrough.R)

## Data visualization

I. Phylogenetic trees (Fig. 1)
* [Phylogenetic trees for mtDNA and NGS data](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/Fig1_phylotrees.R)

II. NGSAdmix analyses (Figs. S1 and S2)
* [Structure-style plots for data (Fig. S1)](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/FigS1_NGSadmix.R)
* [Pie charts plotted on map illustrating NGSadmix proportions (Fig. S2)](https://github.com/eachambers/pantherophis_spp_delim/blob/main/scripts/FigS2_popgen.R)
