# Triatominae_PhyloFields
Functions used in Ceccarelli et al. "Phylogenetic structure of geographic co-occurrence among new world Triatominae species, vectors of chagas disease" Accepted in Journal of Biogeography

File description:

File `ConvexHull.R`: Function to create range map polygons from ocurrence points (from the alphahull R package; Pateiro-Lopez & Rodriguez-Casal 2016). This file should be saved in the working directory, it will be sourced from there (see file Rscripts_PhyloFields_Ceccarellietal2019.R)


File `div_phylo_fields.R`: This file contains the following functions: 

`divfield_pams`: creates and saves (in WD) species' diversity field as a presence-absence matrix of species co-occurring within the a focal species' range.

`sp_assem`: creates and saves focal species' assemblages based on the results (tables) from the `divfield_pams`function.

`psv_modified`: modified version of the function 'psv' from the 'picante' (Kembel et al. 2018 )R package to calculate the phylogenetic structure of species co-ocurring within a focal species' range. 

`phylofieds_psvs`: iterates the 'psv_modified' function for all focal species, thus calculating each focal species' phylogenetic field. 

File `randomrange_pf2.R`: Function to create random phylogenetic fields (randomly generating species assemblages within focal ranges and calculating their phylogenetic structure).

File `compare_PFs.R`: Function to compare observed vs random phylogenetic fields of focal species. 

File `Rscripts_PhyloFields_Ceccarellietal2019.R`: Master script including all analyses. The above functions are to be used within this script. Additional analysis are also included (phylogenetic generalized least squares regressions, using the 'pgls' function from the 'caper' (Orme et al. 2018) R package).
