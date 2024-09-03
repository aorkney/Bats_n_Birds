# Bats_n_Birds

''evolutionary integration of fore- and hindlimb proportions within the bat wing membrane inhibits ecological adaptation compared to birds.''
Orkney, Boerma, Hedrick, (2024)

Code repository for Orkney, Boerma and Hedrick 2024, run in R version 4.2.3 This repository contains commented codes to enable readers to reproduce the analyses conducted in Orkney, Boerma and Hedrick, 2024.

These scripts were run on a machine with the following specifications: Machine specs: 12th Gen Intel(R) Core(TM) i9-12900K 3.20 GHz 128 GB (128 GB usable) 64-bit operating system, x64-based processor

However, no non-standard hardware is required to reproduce the analyses here.

Installation guide: See https://www.r-project.org/ for guidance installing and using R, including expected installation times and the installation of dependant packages.



Data: 

Csize.22.10.2022.RData

This is an array of centroid sizes of landmark constellations, approximating the volumes of several bones
in the avian skeleton for the bird species considered here. 
The source landmark constellations were originally compiled by Bjarnason et al., 2021
(https://doi.org/10.18563/journal.m3.125)
This file will be generated during data preparation and is not supplied in this GitHub repository.

coords.22.10.2022.RData

This object refers to original landmark constellations prepared from Bjarnason et al., 2021. 
(https://doi.org/10.18563/journal.m3.125)
This file will be generated during data preparation and is not supplied in this GitHub repository.

tree.22.10.2022.RData

This is a phylogenetic tree pruned from Prum et al., 2015 (https://doi.org/10.1038/nature15697), 
representing the evolutionary relationships between the bird species considered here.
This file will be generated during data preparation and is not supplied in this GitHub repository. 

tree.names.22.10.2022.RData

This is a vector of species names, which matches the birds in this study to congeners on Prum et al.'s phylogeny. 

masses.22.10.2022.RData

This is a vector of representative body masses for the avian species considered here, 
aggregated from the Cornell handbook of the birds of the world. (https://birdsoftheworld.org/bow/home)

flight_masses_22_10_2022_plus_A.csv

standardised_foot_scores_11_07_2023.csv

These two spreadsheets are standardised presence (1) or absence (0) scores for ecological activites for the 
bird species considered in this study. The flight style scores were expanded from the scheme of Taylor & Thomas
(Taylor, G. & Thomas, A. Evolutionary Biomechanics (Oxford Univ. Press, 2014)), see Orkney et al., 2021 (https://doi.org/10.1038/s41559-021-01509-w)
The foot use scores were originally published in Orkney et al., 2021, and were compiled by Brigit C Tronrud. 

Bird_landmark_info_23Nov2019.csv
This is a metadata file containing information required to extract data from the Bjarnason Landmark deposition.

Eco_meta_data.csv
This is a metadata file containing bird genus names and masses require for data preparation.
(it also contains a variety of other information but it is not needed here)

Bat_CT_process_list_Andrew_only_upload_copy.csv

This is a spreadsheet of metadata related to bat skeletal material that was collected and processed
as part of this study.

Bat_eco_metadata.csv

This is a spreadsheet of presence (1) or absence (0) scores for ecological activites for the bat
species considered in this study. Unknown classifications (?) are treated as absence. 

eco_bibtex.txt

This is a compendium of supporting BibTeX format citations for 'Bat_eco_metadata.csv'.

humerus_array_nov_08_2023.RData

radius_array_nov_08_2023.RData

handwing_array_nov_08_2023.RData

femur_array_nov_08_2023.RData

tibia_array_nov_08_2023.RData

These objects are landmark constellations representing bat skeletal forms, designed and placed by Orkney between 2022 and 2023, 
employed under the Hedrick Lab and leveraging datasets and assets from multiple museums, with extensive support from Boerma.

chiroptera.no_outgroups.absolute.tre 

Bat phylogeny derived from (shi & Rabosky, 2015; https://doi.org/10.1111/evo.12681)
This phylogeny is not provided in this GitHub repository.

Scripts: 

Data preparation:
Bird_data_preparation_08_30_2024.R (Run time under 5 minutes)
This script loads and prepares datasets for subsequent analysis. 
Bird landmark data should be sourced from: https://doi.org/10.18563/journal.m3.125
A family tree describing bird relationships to each other should be sourced from: https://doi.org/10.5281/zenodo.28343
The script will require the following additional files: Bird_landmark_info_23Nov2019.csv ; Eco_meta_data.csv

Fig 1: Figure_1_final.R (Run time on aforementioned machine specifications is 4 minutes)
This figure produces plots that map the pairwise evolutionary covariances between the volumes/centroid sizes (this can be referred to as 'proportions' for brevity) of bones across
bat and bird phylogeny. 
In effect, this figure shows which bones tend to evolve together or independently of one another, in the bird and bat body plan respectively. Body regions that evolve together 
in a consorted way, independently of other body regions, are known as 'evolutionary modules'. 


Fig 2: Pairwise_integration_handwing_resolved.R (Run time on aforementioned machine specifications is 2 minutes)
This figure further probes evolutionary covariances of skeletal proportions within the bat body plan, exploring whether additional insight is gained by subdividing the bat handwing
into its major supporting digits. 

Fig 3: Polar_ecology_plot.R (Run time on aforementioned machine specifications is 11 minutes)
This figure maps the strength of evolutionary covariance between skeletal proportions and multivariate ordinations of flight-style and foot-use/roosting activities across birds and bats, 
in order to explore whether there is a coincidence of evolutionary modules in the bird and bat skeletons with the skeletal regions that exhibit the strongest adaptive responses to distinct
ecological activities. 

Fig 4: Phenogram_pale_core.R (Run time on aforementioned machine specifications is 1 minute)
This figure produces an ancestral state reconstruction of the gross proportions of the appendicular skeleton of birds and bats, and computes indices of the accumulation of phenotypic disparity through time. 
This contributes to discussion, which explores whether differences in evolutionary modules and regionalised ecological adaptation in the bird and bat skeleton may have influenced their gross evolutionary dynamics. 

The aesthetic appearance of plots may differ slightly to the publication versions. 


Supporting materials:

Annotated llustrations of bat humerus, radius, femur and tibia landmark constellations employed to proxy bone sizes:

Humerus_scheme_07_22_2024.png

Radius_scheme_07_22_2024.png

Femur_scheme_07_22_2024.png

Tibia_scheme_07_22_2024.png

Supporting landmark descriptions, with numbers matching annotated images:

Humerus_landmark_descriptions_20_07_2024.txt

Radius_landmark_descriptions_21_07_2024.txt

Femur_landmark_descriptions_21_07_2024.txt

Tibia_landmark_descriptions_21_07_2024.txt




