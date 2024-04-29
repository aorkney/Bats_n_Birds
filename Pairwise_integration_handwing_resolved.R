# The purpose of this script is to
# import the bat dataset and prepare it for analysis...
# We will, along with preparing a block of data that represents 
# overall handwing centroid size, produce a block of data
# that represents:
# 1) Thumb claw length
# 2) Second finger length
# 3) Fourth finger length
# Evolutionary integration will then be computed across pairs of skeletal
# element relative sizes (corrected for body mass), across the bat body plan, 
# including major digit lengths within the handwing.

library(geomorph) # 4.0.5
library(ape) # 5.7-1
library(nlme) # 3.1-162

setwd()
# Set the directory to the location of the landmark data

load('humerus_array_nov_08_2023.RData')
load('radius_array_nov_08_2023.RData')
load('handwing_array_nov_08_2023.RData')
load('femur_array_nov_08_2023.RData')
load('tibia_array_nov_08_2023.RData')
# Load the original bat skeleton landmark arrays. These include unevenly spaced sliding curve landmarks

total.thumb <-matrix(NA, length(dimnames(handwing.array)[[3]]), 1)
total.thumb <- as.vector(total.thumb)
names(total.thumb)<-dimnames(handwing.array)[[3]]
# Make an object to receive thumb length data

total.fourth.finger <- total.second.finger <- total.thumb
# The fourth and second fingers (digits V and III) will have the same dimensionality as thumb length

for( i in 1:length(dimnames(handwing.array)[[3]]
)){
	total.thumb[i] <- handwing.array[,,i][4,2]
	total.second.finger[i] <- handwing.array[,,i][12,2]
	total.fourth.finger[i] <- handwing.array[,,i][18,2]
}
# Populate the digit length objects

bat.landmarks<-list()

bat.landmarks[[1]] <- humerus.array
bat.landmarks[[2]] <- radius.array
bat.landmarks[[3]] <- handwing.array
bat.landmarks[[4]] <- femur.array
bat.landmarks[[5]] <- tibia.array
# Compile the bat landmark data into a list

names(bat.landmarks) <- elements <- c('humerus','radius','handwing','femur','tibia')
# name the list elements

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

# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function


GPA.coords <- lapply( GPA.results , get.item , item = "coords" )
GPA.Csize <- lapply( GPA.results , get.item , item = "Csize" )
# Coordinates and centroid sizes have been extrated
# Centroid sizes, which proxy skeletal element volumes (equivalent to 'proportions' for long, thin bones), are stored in 'GPA.Csize' and will be the subject of analysis

GPA.Csize[[6]] <- total.thumb
GPA.Csize[[7]] <- total.second.finger
GPA.Csize[[8]] <- total.fourth.finger
# Add the digit length data to the list of landmark data

names(GPA.Csize)[c(6,7,8)] <- c('thumb','second.finger','fourth.finger')
# Name the new elements approrpriately

setwd()
# Set work directory to the location of the phylogeny by Shi and Rabosky 

bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# This phylogeny far exceeds the number of taxa for which we have landmark constellations
# we must therefore prune the phylogeny  

pruned.tree.bats <- keep.tip(bat.tree,taxa)
# Pruned phylogeny produced

metadata <- read.csv('Bat_CT_process_list_Andrew_only.csv')
# Read bat metadata

masses.bats <- metadata$Mass[ match(taxa,metadata$Shi_match) ]
names(masses.bats)<-taxa
# Define a vector of bat body masses

# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# The next steps in the code will entail removing covariance between skeletal element proportions that is induced by size-scaling (allometry). 
# The bats in the dataset range over a factor of 1000 in body size, which means that it is necessary to remove this signal so that evolutionary relationships
# that transcend differences between very small and very large bats can be perceived. 

library(phytools) # 1.5-1
# The phytools package contains functions with allow us to initiate our estimate of phylogenetic structure of skeletal element sizes

get.residual.Csize <- function( array, masses, phylogeny, taxa ){
	species<-taxa
	allometry.Csize <- list()
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] )
	for(i in 1:length(array) ){
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ), species = taxa  )
		x <- df$Csize
		names(x)<-rownames(df)
		lambda <- phylosig(tree= newphy, x = x[newphy$tip], method='lambda')[[1]]
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( lambda, phy=newphy, form= ~ species ), data=df )
	}
	names(allometry.Csize) <- names(array)
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c )
	return(residual_Csize)
}

# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest. 

allo.Csize.bats <- get.residual.Csize( array = GPA.Csize, masses = masses.bats, phylogeny = pruned.tree.bats, taxa=taxa )
# Which taxa do we wish to investigate? 
# In this example, we will investigate all birds. 

elements <- c('humerus','radius','thumb','second.finger','fourth.finger','femur','tibia')
# A vector the skeletal elements that will be studied. 

pairs <- combn(elements,2)
# Define a series of unique pairs

Csize.int.bats <- cbind(expand.grid(elements,elements), rep(NA,length(elements)^2),rep(NA,length(elements)^2) )
colnames(Csize.int.bats)<-c('bone1','bone2','Z','p')
# Prepare an object to receive estimates of evolutionary covariance strength between different skeletal element proportion pairs across the bat body plan. 
# A permutation-based Z-score effect size and p-value will be extracted from a phylogenetic 2-Block Partial Least Squares approach. 

for(i in 1:(length(pairs)/2)){ # For each unique pair of skeletal elements
		rows <- intersect( grep( paste(Csize.int.bats$bone1,Csize.int.bats$bone2), pattern= pairs[1,i] )  , grep( paste(Csize.int.bats$bone1,Csize.int.bats$bone2), pattern= pairs[2,i] ) ) 
		# Find the rows in the results object that match the skeletal pair
		size.pls <- phylo.integration(allo.Csize.bats[pairs[1,i]][[1]][pruned.tree.bats$tip],allo.Csize.bats[pairs[2,i]][[1]][pruned.tree.bats$tip], phy=pruned.tree.bats)
		# Compute 2-Block Partial Least Squares analysis 
		Csize.int.bats[rows,]$Z <-  size.pls$Z/sqrt(length(pruned.tree.bats$tip))
		Csize.int.bats[rows,]$p <- size.pls$P.value
		# Extract summary statistics
}


# The following analyses are not integral to the manuscript but may interest the user seeking to replicate the study:
# The 'compare.pls' function determines whether differences between estimated evolutionary covariances are significant. 

compare.pls(phylo.integration(allo.Csize.bats['second.finger'][[1]][pruned.tree.bats$tip],allo.Csize.bats['fourth.finger'][[1]][pruned.tree.bats$tip], phy=pruned.tree.bats),
phylo.integration(allo.Csize.bats['femur'][[1]][pruned.tree.bats$tip],allo.Csize.bats['thumb'][[1]][pruned.tree.bats$tip], phy=pruned.tree.bats)
)
# Comparing integration between the bat finger proportions, to integration between the bat femur and thumb proportions.

compare.pls(phylo.integration(allo.Csize.bats['second.finger'][[1]][pruned.tree.bats$tip],allo.Csize.bats['thumb'][[1]][pruned.tree.bats$tip], phy=pruned.tree.bats),
phylo.integration(allo.Csize.bats['femur'][[1]][pruned.tree.bats$tip],allo.Csize.bats['thumb'][[1]][pruned.tree.bats$tip], phy=pruned.tree.bats)
)
# Compairng integration between the bat thumb and finger proportions, to integration between the bat femur and thumb proportions. 

# The thumb is significantly less integrated to the femur than it is to the second finger, and the fingers are significantly more integrated to each other
# than the thumb is to the femur. 


# The results may now be plotted:

library(ggplot2) # 3.4.2
sizedf.bats <- as.data.frame(Csize.int.bats)
sizedf.bats$Z[which(sizedf.bats$p >0.05)]<-0
# Compile data frame for plotting and mask non-significant evolutionary covariance results. 

# Define a colour blind friendly palette for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Done

sizedf.bats[,1]<-as.character(sizedf.bats[,1])
sizedf.bats[,2]<-as.character(sizedf.bats[,2])

sizedf.bats[,1][which(sizedf.bats[,1]=='second.finger')] <- 'digit III'
sizedf.bats[,2][which(sizedf.bats[,2]=='second.finger')] <- 'digit III'

sizedf.bats[,1][which(sizedf.bats[,1]=='fourth.finger')] <- 'digit V'
sizedf.bats[,2][which(sizedf.bats[,2]=='fourth.finger')] <- 'digit V'
# Housekeeping

# This will produce a pairwise covariance plot 'Battenberg plot' depicting the strength of evolutionary covariances between pairs of skeletal proportions across the bat appendicular skeleton. 
plota<-
ggplot(data=sizedf.bats)+
coord_equal(ratio = 1)+
geom_tile(aes(x=bone1, y=bone2, fill=Z))+
scale_y_discrete(limits=sizedf.bats$bone2[rev(c(29,22,15,8,1,36,43))])+
scale_x_discrete(limits=sizedf.bats$bone2[c(29,22,15,8,1,36,43)])+
scale_fill_gradient2(high='red',low='white',midpoint=0.0,na.value='gray',limits=c(0,1))+
theme(axis.text.x=element_text(size=15, angle=90, vjust = 0.5, hjust=1,colour='black',face = "bold"),
                     axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.title.y=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(size=25,hjust=.5), legend.position="right",panel.background = element_blank())+
	labs(fill = expression(paste( italic('Z/'),sqrt(n) ) ))+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=4.5,ymin=0.5,ymax=1.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=3.5,xmax=4.5,ymin=0.5,ymax=1.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=4.5,xmax=5.5,ymin=1.5,ymax=2.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=4.5,xmax=5.5,ymin=1.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=6.5,xmax=7.5,ymin=3.5,ymax=4.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=6.5,xmax=7.5,ymin=3.5,ymax=4.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=5.5,xmax=6.5,ymin=2.5,ymax=3.5,linewidth=3)+
annotate('rect',colour='white',fill=NA,xmin=5.5,xmax=6.5,ymin=2.5,ymax=3.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=5.5,xmax=7.5,ymin=0.5,ymax=2.5,linewidth=3)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=5.5,xmax=7.5,ymin=0.5,ymax=2.5,linewidth=2)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=5.5,ymin=2.5,ymax=7.5,linewidth=3)+
annotate('rect',colour=cbbPalette[6],fill=NA,xmin=0.5,xmax=5.5,ymin=2.5,ymax=7.5,linewidth=2)


# The following lines will re-compute the above analysis within a re-sampled framework, to produce an impression of the range of confidence intervals. 
# Results will be plotted as a function of the gross distance between different pairs of skeletal elements. 

sizedf.bats <- as.data.frame(Csize.int.bats)
# Prepare data object to receive results. 

results.Csize<-matrix(0,1,1)
n<-100 # Number of replicates
	i<-1 
k<-100 # Number of taxa within any single sub-sample

while(dim(results.Csize)[1] < n*length(pairs[1,]) ){


	taxa <- sample(size=k, dimnames(bat.landmarks[['humerus']])[[3]] , replace=F)
	# sub-sample of bat diversity 

	tryCatch( expr= {result <- get.residual.Csize( array = GPA.Csize, masses = masses.bats, phylogeny = keep.tip(pruned.tree.bats,taxa), taxa=taxa )}, error=function(err) NA)
	# Compute the allometric residuals 

	if(exists('result')){

	pairs <- combn(elements,2)
	# Define a series of unique pairs

	Csize.int<-cbind(pairs[1,],pairs[2,],rep(NA,length(pairs[1,])),rep(NA,length(pairs[1,])))
	Csize.int<-data.frame(Csize.int)
	colnames(Csize.int)[c(3,4)]<-c('Z','p')

	for(j in 1:(length(pairs)/2)){
			# False convergence errors can occur, so we need tow rap this in an error condition statement 
			size.pls <- tryCatch( expr= {phylo.integration(result[pairs[1,j]][[1]][keep.tip(pruned.tree.bats,taxa)$tip],result[pairs[2,j]][[1]][keep.tip(pruned.tree.bats,taxa)$tip], phy=keep.tip(pruned.tree.bats,taxa), print.progress=F)},error=function(err) NA)
			# Compute the evolutionary covariance between the selected pair of skeletal element residual sizes, under a 2-Block Partial Least Squares approach
			Csize.int[j,]$Z <-  size.pls$Z[[1]]/(k^.5)
			Csize.int[j,]$p <- size.pls$P.value[[1]]
			# Extract summary statistics
	}

	if(i==1){
		results.Csize<-Csize.int
		i<-2
	} else {
		results.Csize<-rbind(results.Csize,Csize.int)
	}
	
	rm(Csize.int)
	rm(result)
	}
	#Sundry housekeeping steps

	# This while loop will run until the number of desired replicates has completed. 
	print( dim(results.Csize)[1]/( n*length(pairs[1,])))
	# Update the user on progress.
	
}



pair.identity <- paste(results.Csize$X1,results.Csize$X2)
distance.identity <- rep(NA,length(pair.identity))
distance.identity[which(pair.identity=='humerus radius' | pair.identity=='radius handwing' | pair.identity=='femur tibia' | pair.identity=='second.finger fourth.finger' | pair.identity=='thumb second.finger' | pair.identity=='thumb fourth.finger' | pair.identity=='radius fourth.finger'  | pair.identity=='radius thumb'  | pair.identity=='radius second.finger')] <-1
distance.identity[which(pair.identity=='humerus handwing' | pair.identity=='humerus femur' | pair.identity=='humerus thumb' | pair.identity=='humerus second.finger' | pair.identity=='humerus fourth.finger') ] <-2
distance.identity[which(pair.identity=='humerus tibia' | pair.identity=='radius femur' )] <-3
distance.identity[which(pair.identity=='radius tibia' | pair.identity=='handwing femur' | pair.identity=='thumb femur' | pair.identity=='second.finger femur' | pair.identity=='fourth.finger femur') ] <-4
distance.identity[which(pair.identity=='handwing tibia' | pair.identity=='thumb tibia' | pair.identity=='second.finger tibia' | pair.identity=='fourth.finger tibia') ] <-5
# Define the distances between pairs of skeletal elements. 

module.identity <- rep('other',length(pair.identity))
module.identity[grep('thumb',pair.identity)] <- 'thumb'
# Define which pairs of skeletal elements include the thumb. 

col.pal <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Define a colour palette for plotting. 

df <- data.frame(cbind(pair.identity,distance.identity,results.Csize))
df$Z<-as.numeric(df$Z)
df$p<-as.numeric(df$p)
topo.bats <- ggplot(data=df,aes(x=distance.identity,y=Z, fill=module.identity, shape=module.identity))+
geom_jitter(size=3,alpha=0.5)+
scale_fill_manual(values=col.pal[c(1,4)])+
scale_shape_manual(values=c(21,24))+
lims(y=c(-0.2,1))+
labs(x='topological distance',y=expression(paste( italic('Z/'),sqrt(n) ) ), shape='', fill='' )+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour='black',face = "bold"),
	axis.title.y.left=element_text(size=15,colour='black',face = "bold"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	plot.title=element_text(size=25,hjust=.5),
	legend.position='bottom',
				legend.title=element_text(size=15),
				legend.text=element_text(size=15))
# Produce a scatter plot of the results. 

library(ggpubr) # 0.6.0
# Package for arranging subplots
					     
dev.new(width=25,height=15,unit='cm')
# Make new image
ggarrange(plota, topo.bats, ncol=2,labels=c('a','b'), font.label = list(size = 30, color = "black", face = "bold"), hjust=c(-3,3), widths=c(1,0.8))
# Arrange subplots
setwd()
# Specify directory to save figure
ggsave(filename='Figure_3_03_26_2024.pdf', width = 30,
  height = 18,
  units = c( "cm"),
  dpi = 300)
# Save figure

# Script concludes
