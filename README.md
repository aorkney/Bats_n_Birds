# Bats_n_Birds

''Evolutionary integration of fore- and 
hindlimb proportions within the bat wing membrane 
inhibits ecological adaptation compared to birds.''
Orkney, Boerma, Hedrick, (2024)

Data: 

Csize.22.10.2022.RData
This is an array of centroid sizes of landmark constellations, approximating the volumes of several bones
in the avian skeleton for the bird species considered here. 
The source landmark constellations were originally compiled by Bjarnason et al., 2021
(https://doi.org/10.18563/journal.m3.125)

tree.22.10.2022.RData
This is a phylogenetic tree pruned from Prum et al., 2015 (https://doi.org/10.1038/nature15697), 
representing the evolutionary relationships between the bird species considered here.

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

Bat_CT_process_list_Andrew_only.csv
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
These objects are 

Fig 1: Pairwise_integration_across_bodyplans.R

Fig 2: Pairwise_integration_handwing_resolved.R

Fig 3: Polar_ecology_plot.R

Fig 4: Phenogram_pale_core.R
