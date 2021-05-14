# pantherophis_spp_delim

The following is code to generate figures from [Marshall et al. (2021)](https://www.sciencedirect.com/science/article/abs/pii/S1055790321001275). Raw data files are provided on Dryad [here](LINKHERE).

## Data visualization

I. Phylogenetic trees (Fig. 1)
* [Phylogenetic trees for mtDNA and NGS data](https://github.com/eachambers/pantherophis_spp_delim/blob/main/Fig1_phylotrees.R)

II. NGSAdmix analyses (Figs. S1 and S2)
* [Structure-style plots for data (Fig. S1)](https://github.com/eachambers/pantherophis_spp_delim/blob/main/FigS1_NGSadmix.R)
* [Pie charts plotted on map illustrating NGSadmix proportions (Fig. S2)](https://github.com/eachambers/pantherophis_spp_delim/blob/main/FigS2_popgen.R)

III. Data analysis files
* [Split iPyrad .loci files and save separately](https://github.com/eachambers/pantherophis_spp_delim/blob/main/Split_loci_file.ipynb)
* [Generate .phylip file for each locus for input into RAxML-ng](https://github.com/eachambers/pantherophis_spp_delim/blob/main/processing_loci_files.R)
* [Count number of individuals per locus for input into PHRAPL](https://github.com/eachambers/pantherophis_spp_delim/blob/main/count_pop_loci.R)
* [Subsample pre-loaded PHRAPL MigrationArray](https://github.com/eachambers/pantherophis_spp_delim/blob/main/subsample_migration_array.R)
* [Additional data processing and generate PHRAPL input file](https://github.com/eachambers/pantherophis_spp_delim/blob/main/PHRAPL_workthrough.R)
