# The purpose of this script is to produce pairwise evolutionary covariance/
# 'Battenberg' plots and topological distance plots that portray
# modular schemes of skeletal evolutionary organisation in birds and bats


# Set the work directory
setwd()
# The user will need to customise this

library( geomorph ) # 4.0.5
# Package for macroevolutionary analyses
library( ape ) # 5.7-1
# Package for managing phylogenetic trees
library( nlme ) # 3.1-162
# Package for computing phylogenetic regressions

# Load some requisite datasets
load('tree.22.10.2022.RData') # Avian tree pruned from Prum et al., 2015
load('tree.names.22.10.2022.RData') # Congener names 
load('masses.22.10.2022.RData') # Bird body masses aggregated from Cornell 
load('Csize.22.10.2022.RData') # Centroid sizes of Bjarnason et al., 2021's birds
# Done 

# Ensure names match the phylogeny. 
names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done

Csize.birds<-GPA.Csize

# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

library(phytools) # 1.5-1
# Package with functions that compute phylogenetic signal indices

get.residual.Csize <- function( array, masses, phylogeny, taxa ){
	species<-taxa
	allometry.Csize <- list()
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] )
	for(i in 1:length(array) ){
		lambda <- phylosig(tree= newphy, x=log10( array[[ i ]][taxa]),  method='lambda')[[1]]
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ), species = taxa  )
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( lambda, phy=newphy, form= ~ species ), data=df )
	}
	names(allometry.Csize) <- names(array)
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c )
	return(residual_Csize)
}
# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest. 


allo.Csize.birds <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=names(masses) )
# Compute allometric residuals of centroid sizes. 

names(allo.Csize.birds)[9] <- 'handwing'
names(allo.Csize.birds)[12] <- 'tibia'
# Set the names of the carpometacarpus and tibiotarsus to match their mammalian homologues 

masses.birds <- masses
pruned.tree.birds <- pruned.tree

# The bird dataset has been loaded. 

load('humerus_array_nov_08_2023.RData')
load('radius_array_nov_08_2023.RData')
load('handwing_array_nov_08_2023.RData')
load('femur_array_nov_08_2023.RData')
load('tibia_array_nov_08_2023.RData')
# Load the original landmark arrays describing bat skeletal shape, 
# compiled from scans produced by Boerma, Orkney & Hedrick, 
# and landmarked by Orkney. 
# These constellations include unevenly spaced sliding curve landmarks

bat.landmarks<-list()
bat.landmarks[[1]] <- humerus.array
bat.landmarks[[2]] <- radius.array
bat.landmarks[[3]] <- handwing.array
bat.landmarks[[4]] <- femur.array
bat.landmarks[[5]] <- tibia.array
# Compile constellations into list

names(bat.landmarks) <- elements <- c('humerus','radius','handwing','femur','tibia')
# Compile the landmark data into a list

humerus.sliders <- rbind( define.sliders(c(3, 20:24,9), write.file=F),
define.sliders( c(6,25:31, 8), write.file=F),
define.sliders( c(13,32:41,15), write.file=F),
define.sliders( c(14,42:51, 15), write.file=F),
define.sliders( c(10,52:57, 11), write.file=F),
define.sliders( c(18,58:63, 19), write.file=F) )

radius.sliders <- rbind( define.sliders( c(4,19:23,2) ,write.file=F), 
define.sliders( c(4,24:26,3) ,write.file=F),
define.sliders( c(3,27:30,5) ,write.file=F),
define.sliders( c(5,31:45,2) ,write.file=F),
define.sliders( c(2,36:40,6) ,write.file=F),
define.sliders( c(6,41:44,4) ,write.file=F),
define.sliders( c(5,45:53,7) ,write.file=F),
define.sliders( c(3,54:63,15) ,write.file=F) )

femur.sliders<- rbind( define.sliders(c(2, 17:25,8), write.file=F),
define.sliders( c(6,26:30, 5), write.file=F),
define.sliders( c(10,31:39,16), write.file=F),
define.sliders( c(14,40:48, 15), write.file=F),
define.sliders( c(11,49:53, 12), write.file=F),
define.sliders( c(13,53:58, 12), write.file=F),
define.sliders( c(14,59:68, 1), write.file=F) )

tibia.sliders<- rbind( define.sliders(c(1, 5:9,2), write.file=F),
define.sliders( c(2,10:14, 3), write.file=F),
define.sliders( c(3,15:19,1), write.file=F),
define.sliders( c(1,20:29, 4), write.file=F) )

sliders<-list()
sliders[[1]]<-humerus.sliders
sliders[[2]]<-radius.sliders
sliders[[3]]<-femur.sliders
sliders[[4]]<-tibia.sliders
names(sliders) <- elements[-3]

# Sliding landmark indices have been defined. Note that the handwing does not include sliding landmarks

get.gpa.taxa <- function(landmark.list,sliders,taxa,elements){
	GPA.results <- list() # Define a dumby list to receive aligned landmark constellations
		for( i in 1:length( elements ) )	{ # For each skeletal element
			element <- elements [ i ] # Define the element of interest
			if(length(sliders[[element]])>0){
				GPA.results[[ i ]] <- gpagen( landmark.list[[ element ]][ , , taxa ] , curves=sliders[[element]], approxBE=T ) # Perform Generalised Procrustes Alignment
			} else {
				GPA.results[[ i ]] <- gpagen( landmark.list[[ element ]][ , , taxa ] , approxBE=T ) # Perform Generalised Procrustes Alignment
			}
		}
	names( GPA.results ) <- elements
	return( GPA.results )
}
# This is a function that performs a Generalised Procrustes Alignment on the raw data. 
# The user can specify taxonomic subsets; they may wish to if they suspect some taxa cause a 'pinnochio' effect, for example. 

taxa <- dimnames(bat.landmarks[['humerus']])[[3]]
# For now, let's use all the bat species available

GPA.results <- get.gpa.taxa(bat.landmarks,sliders,taxa,elements)
# The alignemnt has been performed

for(i in 1:length(GPA.results)){
		if(dim(GPA.results[[i]]$coords)[2] <3){ # If a shape constellation is 2D, coerce 
			shape.temp<-GPA.results[[i]]$coords # Define a dumby variable that receives the original shape data for that element
			shape.new<-array(NA, c(dim(shape.temp)[1],3,dim(shape.temp)[3]))
			for(j in 1:(dim(shape.temp)[3])){
				shape.new[,,j] <- cbind(GPA.results[[i]]$coords[,,j], rep(0,dim(GPA.results[[i]]$coords)[1]) )
			}
			dimnames(shape.new)[[1]]<-dimnames(GPA.results[[i]]$coords)[[1]] ; dimnames(shape.new)[[3]]<-dimnames(GPA.results[[i]]$coords)[[3]]
			dimnames(shape.new)[[2]]<-c('X','Y','Z')
			GPA.results[[i]]$coords<-shape.new
			GPA.results[[i]]$consensus <- cbind(GPA.results[[i]]$consensus,rep(0,dim(GPA.results[[i]]$consensus)[1]))
			dimnames(GPA.results[[i]]$consensus)[[2]] <- c('X','Y','Z')
		}
}

GPA.coords <- lapply( GPA.results , get.item , item = "coords" )
GPA.Csize <- lapply( GPA.results , get.item , item = "Csize" )
# Coordinates and centroid sizes have been extracted


bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# Phylogeny sourced from Shi & Rabosky, 2015 for bat relationships
# This phylogeny far exceeds the number of taxa for which we have landmark constellations
# we must therefore prune the phylogeny  
pruned.tree.bats <- keep.tip(bat.tree,taxa)

metadata <- read.csv('Bat_CT_process_list_Andrew_only.csv')
# Metadata file of bat project information

masses.bats <- metadata$Mass[ match(taxa,metadata$Shi_match) ]
names(masses.bats) <- taxa
# A vector of body masses, assembled from Nowak & Walker 1999, and the Encyclopedia of life, has been extracted. 

allo.Csize.bats <- get.residual.Csize( array = GPA.Csize, masses = masses.bats, phylogeny = pruned.tree.bats, taxa=names(masses.bats) )
# Compute the allometric residuals of bat skeletal element centroid sizes. 


pairs <- combn(elements,2)
# Define a series of unique pairs of skeletal elements. 

# We will now explore the strength and significance of evolutionary integration between the allometric residuals of 
# bird and bat skeletal proportions over their respective phylogenetic trees. 

Csize.int.birds <- cbind(expand.grid(elements,elements), rep(NA,length(elements)^2),rep(NA,length(elements)^2) )
colnames(Csize.int.birds)<-c('bone1','bone2','Z','p')

Csize.int.bats <- cbind(expand.grid(elements,elements), rep(NA,length(elements)^2),rep(NA,length(elements)^2) )
colnames(Csize.int.bats)<-c('bone1','bone2','Z','p')

diff.int.animals <- cbind(expand.grid(elements,elements), rep(NA,length(elements)^2),rep(NA,length(elements)^2) )
colnames(diff.int.animals)<-c('bone1','bone2','Z','p')

# Data objects have been defined to receive the results. 
# Integration strength will be computed with a Z-test effect size, comparing the real trait covariance under a 2-Block Partial Least Squares decomposition against
# permuted trait covariance under a Brownian walk. 

bat.models <- list()
bird.models <- list()
for(i in 1:(length(pairs)/2)){
		rows <- intersect( grep( paste(Csize.int.bats$bone1,Csize.int.bats$bone2), pattern= pairs[1,i] )  , grep( paste(Csize.int.bats$bone1,Csize.int.bats$bone2), pattern= pairs[2,i] ) ) 
		# Find the relevant row of the data object to fill
		bat.models[[i]] <- phylo.integration(allo.Csize.bats[pairs[1,i]][[1]][pruned.tree.bats$tip],allo.Csize.bats[pairs[2,i]][[1]][pruned.tree.bats$tip], phy=pruned.tree.bats)
		# Compute a partial least squared decomposition and Z-score
		Csize.int.bats[rows,]$Z <-  bat.models[[i]]$Z/sqrt(length(pruned.tree.bats$tip))
		Csize.int.bats[rows,]$p <- bat.models[[i]]$P.value
		# Assign values to data object

		# Repeat for birds:
		bird.models[[i]] <- phylo.integration(allo.Csize.birds[pairs[1,i]][[1]][pruned.tree.birds$tip],allo.Csize.birds[pairs[2,i]][[1]][pruned.tree.birds$tip], phy=pruned.tree.birds)
		Csize.int.birds[rows,]$Z <-  bird.models[[i]]$Z/sqrt(length(pruned.tree.birds$tip))
		Csize.int.birds[rows,]$p <- bird.models[[i]]$P.value

		# Compute a difference-of-integration Z score. (this is not presented in the manuscript, but is of diagnostic value)
		contrast <- compare.pls(bat.models[[i]],bird.models[[i]])
		diff.int.animals[rows,]$Z <- contrast[[3]][1,2]
		diff.int.animals[rows,]$p <- contrast[[4]][1,2]
}



library(ggplot2) # 3.4.2
# Package for plotting

sizedf.bats <- as.data.frame(Csize.int.bats)
sizedf.bats$Z[which(sizedf.bats$p >0.05)]<-0
sizedf.birds <- as.data.frame(Csize.int.birds)
sizedf.birds$Z[which(sizedf.birds$p >0.05)]<-0
# Convert data objects to data frames and mask non-significant results

diffdf <- as.data.frame(diff.int.animals)
diffdf$Z[which(diffdf$p >0.05)]<-0


# Define a colour blind friendly palette for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Done

plota<-
ggplot(data=sizedf.bats)+
coord_equal(ratio = 1)+
geom_tile(aes(x=bone1, y=bone2, fill=Z))+
scale_y_discrete(limits=sizedf.bats$bone2[rev(c(11,10,1,16,25))])+
scale_x_discrete(limits=sizedf.birds$bone2[(c(11,10,1,16,25))], guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='red',low='white',midpoint=0.0,na.value='gray',limits=c(0,1))+
theme(    panel.background = element_rect(fill='transparent'), #transparent panel bg
    	plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
	axis.text.x=element_text(size=15, vjust=0.3,colour='black',face = "bold"),
                     axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.title.y=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(size=25,hjust=.5), legend.position="right")+
	labs(fill = expression(paste( italic('Z/'),sqrt(n) ) ))+
annotate('rect',colour='black',fill=NA,xmin=1.5,xmax=2.5,ymin=0.5,ymax=1.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=1.5,xmax=2.5,ymin=0.5,ymax=1.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=2.5,xmax=3.5,ymin=1.5,ymax=2.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=2.5,xmax=3.5,ymin=1.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=4.5,xmax=5.5,ymin=3.5,ymax=4.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=4.5,xmax=5.5,ymin=3.5,ymax=4.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=4.5,ymin=2.5,ymax=3.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=3.5,xmax=4.5,ymin=2.5,ymax=3.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=5.5,ymin=0.5,ymax=2.5,linewidth=3)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=5.5,ymin=0.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=2.5,ymax=5.5,linewidth=3)+
annotate('rect',colour=cbbPalette[6],fill=NA,xmin=0.5,xmax=3.5,ymin=2.5,ymax=5.5,linewidth=2)+
ggtitle('Bats')
# Pairwise integration plot across the bat body plan

leg <- cowplot::get_legend(plota)
# Extract legend (you will need to have cowplot installed)


plota<-
ggplot(data=sizedf.bats)+
coord_equal(ratio = 1)+
geom_tile(aes(x=bone1, y=bone2, fill=Z))+
scale_y_discrete(limits=sizedf.bats$bone2[rev(c(11,10,1,16,25))])+
scale_x_discrete(limits=sizedf.birds$bone2[(c(11,10,1,16,25))], guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='red',low='white',midpoint=0.0,na.value='gray',limits=c(0,1))+
theme(    panel.background = element_rect(fill='transparent'), #transparent panel bg
    	plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
	axis.text.x=element_blank(),
                     axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.ticks.x=element_blank(),
				axis.title.y=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(size=25,hjust=.5), legend.position="none")+
	labs(fill = expression(paste( italic('Z/'),sqrt(n) ) ))+
annotate('rect',colour='black',fill=NA,xmin=1.5,xmax=2.5,ymin=0.5,ymax=1.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=1.5,xmax=2.5,ymin=0.5,ymax=1.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=2.5,xmax=3.5,ymin=1.5,ymax=2.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=2.5,xmax=3.5,ymin=1.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=4.5,xmax=5.5,ymin=3.5,ymax=4.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=4.5,xmax=5.5,ymin=3.5,ymax=4.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=4.5,ymin=2.5,ymax=3.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=3.5,xmax=4.5,ymin=2.5,ymax=3.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=5.5,ymin=0.5,ymax=2.5,linewidth=3)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=5.5,ymin=0.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=2.5,ymax=5.5,linewidth=3)+
annotate('rect',colour=cbbPalette[6],fill=NA,xmin=0.5,xmax=3.5,ymin=2.5,ymax=5.5,linewidth=2)+
ggtitle('Bats')
# Regenerate above plot sans legend. 


sizedf.birds2<-sizedf.birds
sizedf.birds2$Z[which(sizedf.birds$Z=='0')]<-NA
# Mask non significant results in birds

plotb<-
ggplot(data=sizedf.birds2)+
coord_equal(ratio = 1)+
geom_tile(aes(x=bone1, y=bone2, fill=Z))+
scale_y_discrete(limits=sizedf.birds$bone2[rev(c(11,10,1,16,25))])+
scale_x_discrete(limits=sizedf.birds$bone2[(c(11,10,1,16,25))], guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='red',low='white',midpoint=0.0,na.value=NA,limits=c(0,1))+
geom_tile(data=sizedf.birds[which(is.na(sizedf.birds$Z)==T),], aes(x=bone1, y=bone2), fill='gray')+
scale_fill_gradient2(high='red',low='white',midpoint=0.0,na.value=NA,limits=c(0,1))+
theme(    panel.background = element_rect(fill='transparent'), #transparent panel bg
    	plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
	axis.text.x=element_text(size=15, vjust=0.3,colour='black',face = "bold"),
                     axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.title.y=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(size=25,hjust=.5), legend.position="none")+
	labs(fill = expression(paste( italic('Z/'),sqrt(n) ) ))+
annotate('rect',colour='black',fill=NA,xmin=1.5,xmax=2.5,ymin=0.5,ymax=1.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=1.5,xmax=2.5,ymin=0.5,ymax=1.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=2.5,xmax=3.5,ymin=1.5,ymax=2.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=2.5,xmax=3.5,ymin=1.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=4.5,xmax=5.5,ymin=3.5,ymax=4.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=4.5,xmax=5.5,ymin=3.5,ymax=4.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=4.5,ymin=2.5,ymax=3.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=3.5,xmax=4.5,ymin=2.5,ymax=3.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=5.5,ymin=0.5,ymax=2.5,linewidth=3)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=5.5,ymin=0.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=2.5,ymax=5.5,linewidth=3)+
annotate('rect',colour=cbbPalette[6],fill=NA,xmin=0.5,xmax=3.5,ymin=2.5,ymax=5.5,linewidth=2)+
ggtitle('Birds')
# Pairwise plot of evolutionary integration across the avian body plan



# The following loop will repeat the above analyses in a sub-sampled framework, and re-plot them in an alternative visualisation, organised by the distance between pairs of bones:

elements <- c('humerus','radius','handwing','femur','tibia')
pairs <- combn(elements,2)
results.Csize.bats<-matrix(0,1,1)
results.Csize.birds<-matrix(0,1,1)
n<-100
	i<-1 
k<-100
while(dim(results.Csize.bats)[1] < n*length(pairs[1,]) ){


	taxa <- sample(size=k, dimnames(bat.landmarks[['humerus']])[[3]] , replace=F)
	# sub-sample of bat diversity 

	tryCatch( expr= {result <- get.residual.Csize( array = GPA.Csize, masses = masses.bats, phylogeny = keep.tip(pruned.tree.bats,taxa), taxa=taxa )}, error=function(err) NA)
	

	if(exists('result')){

	pairs <- combn(elements,2)
	# Define a series of unique pairs


	Csize.int<-cbind(pairs[1,],pairs[2,],rep(NA,length(pairs[1,])),rep(NA,length(pairs[1,])))
	Csize.int<-data.frame(Csize.int)
	colnames(Csize.int)[c(3,4)]<-c('Z','p')

	for(j in 1:(length(pairs)/2)){
			# False convergence errors can occur, so we need tow rap this in an error condition statement 
			size.pls <- tryCatch( expr= {phylo.integration(result[pairs[1,j]][[1]][keep.tip(pruned.tree.bats,taxa)$tip],result[pairs[2,j]][[1]][keep.tip(pruned.tree.bats,taxa)$tip], phy=keep.tip(pruned.tree.bats,taxa), print.progress=F)},error=function(err) NA)
			Csize.int[j,]$Z <-  size.pls$Z[[1]]/(k^.5)
			Csize.int[j,]$p <- size.pls$P.value[[1]]
	}

	if(i==1){
		results.Csize.bats<-Csize.int
		i<-2
	} else {
		results.Csize.bats<-rbind(results.Csize.bats,Csize.int)
	}
	
	rm(Csize.int)
	rm(result)
	}

	print( dim(results.Csize.bats)[1]/( n*length(pairs[1,])))
	
}





i<-1

while(dim(results.Csize.birds)[1] < n*length(pairs[1,]) ){


	taxa <- sample(size=k, names(Csize.birds[['humerus']]) , replace=F)
	# sub-sample of bird diversity 

	tryCatch( expr= {result <- get.residual.Csize( array = Csize.birds, masses = masses.birds, phylogeny = keep.tip(pruned.tree.birds,taxa), taxa=taxa )}, error=function(err) NA)
	

	if(exists('result')){

	pairs <- combn(elements,2)
	# Define a series of unique pairs

	names(result)[9]<-'handwing'
	names(result)[12]<-'tibia'


	Csize.int<-cbind(pairs[1,],pairs[2,],rep(NA,length(pairs[1,])),rep(NA,length(pairs[1,])))
	Csize.int<-data.frame(Csize.int)
	colnames(Csize.int)[c(3,4)]<-c('Z','p')

	for(j in 1:(length(pairs)/2)){
			# False convergence errors can occur, so we need tow rap this in an error condition statement 
			size.pls <- tryCatch( expr= {phylo.integration(result[pairs[1,j]][[1]][keep.tip(pruned.tree.birds,taxa)$tip],result[pairs[2,j]][[1]][keep.tip(pruned.tree.birds,taxa)$tip], phy=keep.tip(pruned.tree.birds,taxa), print.progress=F)},error=function(err) NA)
			Csize.int[j,]$Z <-  size.pls$Z[[1]]/(k^.5)
			Csize.int[j,]$p <- size.pls$P.value[[1]]
	}

	if(i==1){
		results.Csize.birds<-Csize.int
		i<-2
	} else {
		results.Csize.birds<-rbind(results.Csize.birds,Csize.int)
	}
	
	rm(Csize.int)
	rm(result)
	}

	print( dim(results.Csize.birds)[1]/( n*length(pairs[1,])))
	
}


pair.identity <- paste(results.Csize.bats$X1,results.Csize.bats$X2)
distance.identity <- rep(NA,length(pair.identity))
distance.identity[which(pair.identity=='humerus radius' | pair.identity=='radius handwing' | pair.identity=='femur tibia')] <-1
distance.identity[which(pair.identity=='humerus handwing' | pair.identity=='humerus femur') ] <-2
distance.identity[which(pair.identity=='humerus tibia' | pair.identity=='radius femur') ] <-3
distance.identity[which(pair.identity=='radius tibia' | pair.identity=='handwing femur') ] <-4
distance.identity[which(pair.identity=='handwing tibia') ] <-5
module.identity <- rep('cross',length(pair.identity))
module.identity[which(pair.identity== 'humerus radius' | pair.identity=='radius handwing' | pair.identity=='humerus handwing')] <- 'wing'
module.identity[which(pair.identity=='femur tibia')] <- 'leg'
homo.identity <- rep('no',length(pair.identity))
homo.identity[which(pair.identity=='radius tibia' | pair.identity=='humerus femur'  )]<-'yes'

module.identity[which(homo.identity=='yes')] <- 'serial'

col.pal <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


df <- data.frame(cbind(pair.identity,distance.identity,results.Csize.bats))
df$Z<-as.numeric(df$Z)
df$p<-as.numeric(df$p)
topo.bats <- ggplot(data=df,aes(x=distance.identity,y=Z,fill=module.identity, shape=module.identity))+
geom_jitter(size=3,alpha=0.5)+
scale_fill_manual(values=c(col.pal[c(1)],col.pal[c(8)],'white',col.pal[c(3)]) )+
scale_shape_manual(values=c(21:24))+
lims(y=c(-0.2,1))+
ggtitle('Bats')+
scale_x_continuous(position = "top")+
scale_y_continuous(position = "right", limits=c(-0.2,1))+
labs(x='topological distance',y=expression(paste( italic('Z/'),sqrt(n) ) ), shape='module:', fill='module:' )+
  theme(axis.line.x.top=element_line(size = 1, colour = "black"),
	axis.line.y.right=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.top=element_text(size=15,colour='black',face = "bold"),
	axis.title.y.right=element_text(size=15,colour='black',face = "bold"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	plot.title=element_text(size=25,hjust=.5),
	legend.position='top',
				legend.title=element_text(size=15),
				legend.text=element_text(size=15))

mod.leg <- cowplot::get_legend(topo.bats)

topo.bats <- ggplot(data=df,aes(x=distance.identity,y=Z,fill=module.identity, shape=module.identity))+
geom_jitter(size=3,alpha=0.5)+
scale_fill_manual(values=c(col.pal[c(1)],col.pal[c(8)],'white',col.pal[c(3)]) )+
scale_shape_manual(values=c(21:24))+
lims(y=c(-0.2,1))+
ggtitle('Bats')+
scale_x_continuous(position = "top")+
scale_y_continuous(position = "right", limits=c(-0.2,1))+
labs(x='topological distance',y=expression(paste( italic('Z/'),sqrt(n) ) ), shape='module:', fill='module:' )+
  theme(axis.line.x.top=element_line(size = 1, colour = "black"),
	axis.line.y.right=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	#axis.title.x.top=element_text(size=15,colour='black',face = "bold"),
	axis.title.x.top=element_blank(),
	axis.title.y.right=element_text(size=15,colour='black',face = "bold"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	plot.title=element_text(size=25,hjust=.5),
	legend.position='none',
				legend.title=element_text(size=15),
				legend.text=element_text(size=15))



df <- data.frame(cbind(pair.identity,distance.identity,results.Csize.birds))
df$Z<-as.numeric(df$Z)
df$p<-as.numeric(df$p)
topo.birds <- ggplot(data=df,aes(x=distance.identity,y=Z,fill=module.identity, shape=module.identity))+
geom_jitter(size=3,alpha=0.5)+
scale_fill_manual(values=c(col.pal[c(1)],col.pal[c(8)],'white',col.pal[c(3)]) )+
scale_shape_manual(values=c(21:24))+
lims(y=c(-0.2,1))+
ggtitle('Birds')+
scale_y_continuous(position = "right", limits=c(-0.2,1))+
labs(x='topological distance',y=expression(paste( italic('Z/'),sqrt(n) ) ), shape='module:', fill='module:' )+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.right=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour='black',face = "bold"),
	axis.title.y.right=element_text(size=15,colour='black',face = "bold"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	plot.title=element_text(size=25,hjust=.5),
	legend.position='none',
				legend.title=element_text(size=15),
				legend.text=element_text(size=15))


library(ggpubr)
left.pane <- ggarrange(plota,plotb,ncol=1,common.legend=T,legend='left', labels=c('a','c'), font.label = list(size = 30, color = "black", face = "bold"))

right.pane <- ggarrange(topo.bats,topo.birds,ncol=1,common.legend=T,legend='bottom', labels=c('b','d'), font.label = list(size = 30, color = "black", face = "bold"))


library(jpeg)

# img <- readJPEG("Figure_2_11_09_2023c_done_cropped.jpg")


library(cowplot) # version 1.1.1

blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Done

dev.new(width=20,height=15,unit='cm')

plot.with.inset <-
  ggdraw() +
  draw_plot(blank) +
  #draw_image(img,x=0.08,y=0.01,width=0.8,height=.85)+ 
  draw_plot(plota, x = 0.02, y = 0.29, width = .29, height = 1)+
  draw_plot(plotb, x = 0.02, y = -0.285, width = .29, height = 1)+
	draw_plot(leg,x=-0.05,y=0.41,width=0.2,height=0.25)+
  draw_plot(topo.bats, x = 0.67, y = 0.575, width = .32, height = .43)+
  draw_plot(topo.birds, x = 0.67, y = 0.0, width = .32, height = .43)+
	draw_plot(mod.leg,x=0.7,y=0.38,width=0.2,height=0.3)+
geom_text(data=data.frame(cbind(x=c(0.04,0.96,0.04,0.96),y=c(0.98,0.98,0.42,0.42))), aes(x=x,y=y), label=c('a','b','c','d'), size = 30/.pt, color = "black", fontface = "bold")


plot.with.inset <-
  ggdraw() +
  draw_plot(blank) +
  #draw_image(img,x=0.08,y=0.018,width=0.8,height=.9)+ 
  draw_plot(plota, x = 0.02, y = 0.33, width = .29, height = 1)+
  draw_plot(plotb, x = 0.02, y = -0.265, width = .29, height = 1)+
	draw_plot(leg,x=-0.05,y=0.43,width=0.2,height=0.25)+
  draw_plot(topo.bats, x = 0.67, y = 0.575, width = .32, height = .43)+
  draw_plot(topo.birds, x = 0.67, y = 0.0, width = .32, height = .43)+
	draw_plot(mod.leg,x=0.7,y=0.38,width=0.2,height=0.35)+
geom_text(data=data.frame(cbind(x=c(0.04,0.96,0.04,0.96),y=c(0.98,0.98,0.42,0.42))), aes(x=x,y=y), label=c('a','b','c','d'), size = 30/.pt, color = "black", fontface = "bold")


setwd()
ggsave(filename='lead_figure_02_20_2024.pdf', width = 40,
  height = 30,
  units = c( "cm"),
  dpi = 300)