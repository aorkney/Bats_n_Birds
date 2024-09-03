# The purpose of this script is to prepare
# a suite of geometric morphometric landmark constellations that can be used
# to perform the analyses integral to Orkney & Hedrick 2023

# The source data will consist of a dataset of landmark constellations compiled by Bjarnason & Benson 2021, 
# describing the shapes of skeletal elements of 149 bird species, ( https://doi.org/10.18563/journal.m3.125 )
# a phylogeny produced by Prum et al., 2015, (  https://doi.org/10.5281/zenodo.28343 )
# and various metadata:
# Bird_landmark_info_23Nov2019.csv ; Eco_meta_data.csv

# The 'geomorph' package will be needed
# If you do not have this package, run 'install.packages('geomorph') and proceed as guided
require( geomorph )
# Package loaded

# Load the original dataset from Bjarnason and Benson 2021 https://doi.org/10.18563/journal.m3.125
load("D:/Documents/Alex_birds/NEW_analysis/2021_Bjarnason_125_SIdata/Supp 3 Datasets for Bjarnason & Benson 2021 revised submission 26 January 2021/Bird_landmarks_clean_26January2021_Bjarnason&Benson taxon set.RData")
# Data loaded; the file pathway will need to be customised by the user


landmark.info<-read.csv("D:\\Documents\\Alex_birds\\Bird_landmark_info_23Nov2019.csv" , row.names = 1)

# Set the work directory (this may or may not be needed, depending on the location the user has stored data)
setwd("D:/Documents/Alex_birds/")
# Work directory has been set 

# The 'ape' package will be needed
# If you do not have this package, run 'install.packages('ape') and proceed as guided
library(ape)
# Package loaded

# Load the Avian Phylogeny of Prum et al., (2015)
PrumTree <- read.nexus(file = "Avian-TimeTree.tre")
# Phylogeny loaded

# Set the work directory (this may or may not be needed, depending on the location the user has stored data)
setwd('D:/Documents/Alex_birds/NEW_analysis/')
# Work directory has been set 

# Load ancilary metadata (this contains information about bird masses and ecological properties should the user be interested)
metadata <- read.csv( "Eco_meta_data.csv" , row.names = 1 )
# Metadata loaded

# Semilandmark curves are employed in the landmark constellations of Bjarnason & Benson 2021; set their length to minimum
SL.counts <- "min"
# Value assigned

# Define a vector of all skeletal elements described by the landmark constellations of Bjarnason & Benson 2021
elements <- c( "skull" , "mandible" , "scapula", "coracoid", "sternum" , "humerus", "ulna", "radius", "carpometacarpus", "synsacrum" , "femur", "tibiotarsus", "tarsometatarsus" )
# Vector defined

# The following lines will extract the landmark constellations for each bird, for each skeletal element, with semilandmark curves set to minimum length
# Esoteric assignments will be performed to correct spelling mistakes and remove duplicate data
# A Generalised Procrustes Analysis will then be performed to scale, rotate and align the landmark constellations of each species within each bone, 
# so that interspecific difference may be analysed; the commenting style will follow individual lines while this takes place

name_list<-list() # Define a dumby vector
for(i in 1:length(elements)){ # For each skeletal element 
	name_list[[i]]<-dimnames(compiled.bird.landmarks.list[[ elements[ i ] ]][[ SL.counts ]][[ "output.landmarks.array" ]] )[[ 3 ]] # Extract the required dimension names
}
AA<-table(unlist(name_list)) # Compile the dimension names for all skeletal elements into a table
birds <- names( AA )[ AA == length( elements ) ] # Extract the dimension that corresponds to the bird species names
birds <- birds[ birds != "Falco_sparverius_UMMZ_154452" ] # Remove an unwanted name

GPA.results <- list() # Define a dumby list to receive aligned landmark constellations
	for( i in 1:length( elements ) )	{ # For each skeletal element
		element <- elements [ i ] # Define the element of interest
			GPA.results[[ i ]] <- gpagen( compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "output.landmarks.array" ]][ , , birds ] ,  curves = compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "sliders" ]] , ProcD = T ) # Perform Generalised Procrustes Alignment
		}
names( GPA.results ) <- elements # Ensure that the names of the list elements reflect the skeletal elements that they represent

birds[7]<-'Anseranas semipalmata NHMUK 1852.7.22.1' # Correct a spelling mistake in the metadata
birds[10]<-'Archilochus_colubris_FMNH_484738' # Correct a spelling mistake in the metadata
birds[12]<-'Choriotis_kori_UMMZ_210444' # Correct a spelling mistake in the metadata
birds[13]<-'Arenaria_interpres_FMNH_313992' # Correct a spelling mistake in the metadata
birds[30]<-'Pipra_erythrocephala_UMMZ_157210' # Correct a spelling mistake in the metadata
birds[31]<-'Chaetura_brachyrua_UMMZ_157689' # Correct a spelling mistake in the metadata
birds[32]<-'Charadris_vociferus_FMNH_470173' # Correct a spelling mistake in the metadata
birds[37]<-'Larus_novae-hollandiae_UMZC_274.C_' # Correct a spelling mistake in the metadata
birds[44]<-'Columbina_minuta_FMNH_289160' # Correct a spelling mistake in the metadata
birds[48]<-'Chloropipio_holochlora_FMNH_288169' # Correct a spelling mistake in the metadata
birds[66]<-'Podica_senegalensis_UMZC_209A_' # Correct a spelling mistake in the metadata
birds[75]<-'Leptotila_rufazilla_FMNH_318659' # Correct a spelling mistake in the metadata
birds[84]<-'Myiobius_varbatus_FMNH_386812' # Correct a spelling mistake in the metadata
birds[105]<-'Poecile_atricapillis_FMNH_504323' # Correct a spelling mistake in the metadata
birds[107]<-'Probosciger_aterriumus_NHMUK_S2006.15.5' # Correct a spelling mistake in the metadata
birds[117]<-'Rhamphastos_ambiguus_NHMUK_S2002.21' # Correct a spelling mistake in the metadata
birds[124]<-'Rupicola_peruviana_UMMZ_119210' # Correct a spelling mistake in the metadata

# The stated errands have concluced: A Generalised Procrustes Alignment has been performed on each skeletal element; birds absent from the wider meta data have been removed;
# All of the rotated constellations host 149 bird species

# A sundry function will be defined to extract coordinate and centroid size geometric morphometric data from the aligned landmark constellations
get.item <- function( X , item ) { X[[ item ]] }	
# Function has been defined

# The coordinates and centroid sizes will now be extracted
GPA.coords <- lapply( GPA.results , get.item , item = "coords" )
GPA.Csize <- lapply( GPA.results , get.item , item = "Csize" )
# Coordinates and centroid sizes have been extrated

# The names of species and the tip.labels of the Prum et al., 2015 phylogeny need to correspond
tree_genera <- as.character( metadata[ birds , "Prum_genus" ] )
tree_names <- PrumTree$tip.label[ match( tree_genera , sub( "_.*" , "\\1" , PrumTree$tip.label ) ) ]
# Correspondance has been achieved

# Species that are not needed need to be pruned fro the Prum et al., 2015 phylogeny
pruned.tree <- drop.tip( PrumTree , PrumTree$tip.label[ !PrumTree$tip.label %in% tree_names ] )
# The phylogeny has been pruned

# A vector of masses for the study bird species will be produced by averaging reported male and female body masses sourced from the wider literature
masses <- apply( metadata[ birds , c( "Mass_F_.hbw_alive." , "Mass_M_.hbw_alive." ) ] , 1 , mean )
# A vector of body masses has been defined

# The data that have been compiled must now be saved before downstream analyses are performed

# Set the work directory as desired ( the user will need to customise this )
setwd()
# Work directory should be set by user

# Save data files; first argument corresponds to object, second argument to the saved file name
save(masses, file='masses.22.10.2022.RData')
save(pruned.tree,file='tree.22.10.2022.RData')
save(GPA.Csize,file='Csize.22.10.2022.RData')
save(GPA.coords,file='coords.22.10.2022.RData')
save(tree_names,file='tree.names.22.10.2022.RData')
# Files have been saved
