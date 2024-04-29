# The purpose of this script is to produce a violin plot illustrating the strength of
# evolutionary covariance between ecological variety in flight-style and foot-use/roosting and anatomical properties of the 
# bird and bat skeletons. 
# Skeletal relative proportions and various indices of shape or kinematic ratios
# will be investigated as proxies of skeletal anatomy. 
# Ecological variety is represented by multivariate ordinations of presence/absence scores for ecological activities.

# The first step undertaken will be to remove size-dependent scaling in skeletal proportions from the bat skeletal proportion data.
# It is necessary to define functions to achieve this in a phylogenetic context. 
# It is necessary to remvove size-dependent scaling because the study bats range over a factor of 1000 in body mass, which means
# that covariances in skeletal proportions could easily be overwhelmed by size-dependent scaling, obfuscating evolutionary covariances
# between skeletal proportions across bats. 

# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

library(phytools) # 1.5-1
# The phytools package contains functions to produce an initial estimate of phylogenetic structure in the skeletal proportion data
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
# This function takes an array of skeletal element proportions, species body masses, a phylogeny and a vector of species. 
# It will compute a generalised least squares fit model of size-dependent scaling between bones and body mass in a phylogenetic context and return the residuals for downstream analysis. 

library( geomorph ) # 4.0.5
library( ape ) # 5.7-1
library( nlme ) # 3.1-162
library(cluster) # 2.1.4

setwd()
# Set the directory to the location of the landmark data

load('humerus_array_nov_08_2023.RData')
load('radius_array_nov_08_2023.RData')
load('handwing_array_nov_08_2023.RData')
load('femur_array_nov_08_2023.RData')
load('tibia_array_nov_08_2023.RData')
# Load the original landmark arrays. These include unevenly spaced sliding curve landmarks

bat.landmarks<-list()

bat.landmarks[[1]] <- humerus.array
bat.landmarks[[2]] <- radius.array
bat.landmarks[[3]] <- handwing.array
bat.landmarks[[4]] <- femur.array
bat.landmarks[[5]] <- tibia.array

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

# The following lines standardise the original landmark contellation data formats to 3D arrays
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
# Coordinates and centroid sizes have been extracted; centroid sizes are the proxy for skeletal proportions that will be subject to downstream analysis.

bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# This phylogeny far exceeds the number of taxa for which we have landmark constellations
# we must therefore prune the phylogeny  

pruned.tree.bats <- keep.tip(bat.tree,taxa)

metadata <- read.csv('Bat_CT_process_list_Andrew_only.csv')
masses.bats <- metadata$Mass[ match(taxa,metadata$Shi_match) ]
names(masses.bats) <- taxa
# Read bat metadata, extract body masses

metadata <- read.csv('Bat_eco_metadata.csv')
# Read bat ecological data
order <- match(taxa,metadata$Shi) 
# A vector of row identities that match the ecology dataset to the sample bats for downstream analysis
eco <- metadata[order,]
# The ecological score matrix has been ordered and prepared

binary.flight.scores <- eco[,c(21:35)]
# Extract flight scores from ecological data
binary.flight.scores<-as.matrix(binary.flight.scores)
# Prepare the data as a matrix to make it easy to index. 
binary.flight.scores[which(binary.flight.scores=='?' | binary.flight.scores=='' )] <- 0 
binary.flight.scores[grep(binary.flight.scores,pattern='\\?')] <- 0 
# Set unknown scores to absences 

binary.flight.scores <- apply(binary.flight.scores,2,FUN=as.numeric)
# I am making sure the data class is numeric
rownames(binary.flight.scores)<-taxa
# Set the row names of the flight style matrix to our taxa 
binary.flight.scores <- binary.flight.scores[,-c(13,15)]
# prune unnecessary flight style categorisations (too few bats practiced these styles in the sample data)

binary.roost.scores <- eco[,c(6:16)]
# Extract roost score data 
binary.roost.scores<-as.matrix(binary.roost.scores)
# Prepare the data as a matrix to make it easy to index. 
binary.roost.scores[which(binary.roost.scores=='?' | binary.roost.scores=='' | is.na(binary.roost.scores))] <- 0 
binary.roost.scores[grep(binary.roost.scores,pattern='\\?')] <- 0 
# Set unknown scores to absences

binary.roost.scores <- apply(binary.roost.scores,2,FUN=as.numeric)
# I am making sure the data class is numeric
rownames(binary.roost.scores)<-taxa
# Set the row names of the roosting matrix to our taxa 

total.thumb <-matrix(NA, length(taxa), 1)
total.thumb <- as.vector(total.thumb)
names(total.thumb)<-taxa
# Produce a thumb length data object

total.fourth.finger <- total.second.finger <- total.thumb

for( i in 1:length(taxa)){
	total.thumb[i] <- handwing.array[,,i][4,2]
	total.second.finger[i] <- handwing.array[,,i][12,2]
	total.fourth.finger[i] <- handwing.array[,,i][18,2]
}
# Produce data objects for the lenghs of other major digits


anatomy <- c(names(GPA.Csize),'Handwing.3D','gross')
# Prodyce a vector of skeletal element and multi-element data object names for the bat skeleton, the ecological adaptive signal in which we will explore. 

ecology <- c('flight','roost')
# Define a vector of the ecological category types we will explore. 


# First, we are going to compute an analysis with all taxa considered, in order guage the significance
# of relationships. 
# Thereafter, we will produce a confidence interval of Z values (effect size scores) by assessing 
# various subsamples. 

p.total <- list()
for(i in 1:length(ecology)){
	p.total[[i]] <- matrix(NA, 1, length(anatomy))
	colnames(p.total[[i]])<-anatomy
}
names(p.total)<-ecology
# This object receives the p-value of 2-Block Partial Least Square analyses for evolutionary covariance between skeletal element/multi-element data objects and multivariate ordinations of flight or roosting style ecology. 


allo.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses.bats, phylogeny = pruned.tree.bats, taxa=taxa  )
allo.Csize.bats <- allo.Csize
# Compute allometric residuals of appendicular element proportions in bats
gross <- cbind( GPA.Csize$humerus[taxa]+GPA.Csize$radius[taxa]+total.second.finger[taxa],
total.fourth.finger[taxa],total.thumb[taxa])
gross <- gross/rowSums(gross)
# Compute gross wing proportions- a multivariate representation of bat wing shape, comprising total wing length, digit V length and thumb length. 
Handwing.3D <-GPA.coords$handwing[,,taxa]
# Extract the internal proportions of the bat handwing- a multivariate data object representing handwing shape. Note that no allometric correction is performed on the gross wing proportions or handwing shape variables. 
allo.Csize.bats[[6]] <- Handwing.3D
names(allo.Csize.bats)[6] <- c('Handwing.3D') 
allo.Csize.bats[[7]] <- gross
names(allo.Csize.bats)[7] <- c('gross') 
# Sundry steps to add the multivariate data objects to the total landmark data.

mat <- binary.flight.scores[taxa,]
# Define an ordered matrix for the flight style scores
var<-apply(mat,2,sd)
# Check that columns are not invariant 
mat <- mat[,which(var>0)]
# Remove invariant columns 
flight.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
# transform the flight scores into a distance matrix 
if(length(which(is.na(flight.dist)==T))>0){
	flight.dist[is.na(flight.dist)]<-0
}
# adjust the diagonal to 0 
flight.pcoa <- cmdscale(flight.dist, eig=T)
# Find the explained variances of the eigenvectors of this object
k<-max(which(flight.pcoa$eig/sum(flight.pcoa$eig)>=0.05))
# Consider only those axes which explain more than 5% of variance
flight.pcoa <- cmdscale( flight.dist, k=k  ,eig = T)
# Compute a principal coordinate analysis- a multivariate ordination of the flight style data, constrained to the axes of interest

# Repeat the above for roosting variety
mat <- binary.roost.scores[taxa,]
var<-apply(mat,2,sd)
mat <- mat[,which(var>0)]
roost.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
if(length(which(is.na(roost.dist)==T))>0){
	roost.dist[is.na(roost.dist)]<-0
}
roost.pcoa <- cmdscale(roost.dist, eig=T)
k<-max(which(roost.pcoa$eig/sum(roost.pcoa$eig)>=0.05))
roost.pcoa <- cmdscale( roost.dist, k=k  ,eig = T)

pcoas <- list(flight.pcoa,roost.pcoa)
names(pcoas)<-ecology
# Sundry step to home the multivariate ecology objects in one object

# The following loops will compute 2-Block Partial Least Squares analyses of evolutionary covariance strength between the anatomical variables and the ecological variables in all anatomy-ecology pairwise combinations. 
i<-1
for(j in 1: length(anatomy)){
	for(k in 1:length(ecology)){
		if( is.null(dim(allo.Csize.bats[[anatomy[j]]])) ){
			tryCatch( expr= {model <- phylo.integration(allo.Csize.bats[[anatomy[j]]][pruned.tree.bats$tip], pcoas[[ecology[k]]]$points[pruned.tree.bats$tip,], phy=pruned.tree.bats,print.progress=F ) }, error=function(err) NA)
		} else  if( length(dim(allo.Csize.bats[[anatomy[j]]]) )<3 ) {
			tryCatch( expr= {model <- phylo.integration(allo.Csize.bats[[anatomy[j]]][pruned.tree.bats$tip,], pcoas[[ecology[k]]]$points[pruned.tree.bats$tip,], phy=pruned.tree.bats,print.progress=F ) }, error=function(err) NA)			
		} else {
			tryCatch( expr= {model <- phylo.integration(allo.Csize.bats[[anatomy[j]]][,,pruned.tree.bats$tip], pcoas[[ecology[k]]]$points[pruned.tree.bats$tip,], phy=pruned.tree.bats,print.progress=F ) }, error=function(err) NA)
		}
		 p.total[[k]][i,j] <- model$P.value[[1]] #results.Z[[k]][i,j] <- model$Z[[1]] ; # Only the p-value will be extracted, because the purpose of this loop is to assess significance. 
	}
}

# p.total demonstrates that the gross proportions of the bat wing are consistently significantly associated to all ecological categories. 
# Gross wing proportions and handwing proportions are both strongly associated with flight style
# Gross wing proportions, femur and tibia size and the size of the handwing are all significantly associated to roosting style variety. 


# The following loops will perform the analyses in a subsetted re-sampled framework across bats, and Z-scores will be extracted to represent the strength of evolutionary covariances between anatomical and ecology pairs. 
results.Z <- list()
for(i in 1:length(ecology)){
	results.Z[[i]] <- matrix(NA, 100, length(anatomy))
	colnames(results.Z[[i]])<-anatomy
}
names(results.Z)<-ecology

results.p <- results.Z 


for(i in 1: 100){

	a<-1
	while(a<2){
		subsample <- sample(taxa,100)
		
		tryCatch( expr= {allo.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses.bats, phylogeny = pruned.tree.bats, taxa=subsample  ) }, error=function(err) NA)
		if(exists('allo.Csize')){
			a<-2
		}
		sub.tree <- keep.tip(pruned.tree.bats, subsample)
		
	}
	allo.Csize.bats <- allo.Csize

	gross <- cbind( GPA.Csize$humerus[subsample]+GPA.Csize$radius[subsample]+total.second.finger[subsample],
	total.fourth.finger[subsample],total.thumb[subsample])
	gross <- gross/rowSums(gross)

	#silhouette <-GPA.coords$handwing[,,subsample]
	Handwing.3D <-GPA.coords$handwing[,,subsample]

	
	allo.Csize.bats[[6]] <- Handwing.3D
	names(allo.Csize.bats)[6] <- c('Handwing.3D') 

	allo.Csize.bats[[7]] <- gross
	names(allo.Csize.bats)[7] <- c('gross') 

	mat <- binary.flight.scores[subsample,]
	var<-apply(mat,2,sd)
	mat <- mat[,which(var>0)]
	flight.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
	if(length(which(is.na(flight.dist)==T))>0){
		flight.dist[is.na(flight.dist)]<-0
	}
	flight.pcoa <- cmdscale(flight.dist, eig=T)
	k<-max(which(flight.pcoa$eig/sum(flight.pcoa$eig)>=0.05))
	flight.pcoa <- cmdscale( flight.dist, k=k  ,eig = T)


	mat <- binary.roost.scores[subsample,]
	var<-apply(mat,2,sd)
	mat <- mat[,which(var>0)]
	roost.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
	if(length(which(is.na(roost.dist)==T))>0){
		roost.dist[is.na(roost.dist)]<-0
	}
	roost.pcoa <- cmdscale(roost.dist, eig=T)
	k<-max(which(roost.pcoa$eig/sum(roost.pcoa$eig)>=0.05))
	roost.pcoa <- cmdscale( roost.dist, k=k  ,eig = T)

	pcoas <- list(flight.pcoa,roost.pcoa)
	names(pcoas)<-ecology

	for(j in 1: length(anatomy)){
		for(k in 1:length(ecology)){
			if( is.null(dim(allo.Csize.bats[[anatomy[j]]])) ){
				tryCatch( expr= {model <- phylo.integration(allo.Csize.bats[[anatomy[j]]][sub.tree$tip], pcoas[[ecology[k]]]$points[sub.tree$tip,], phy=sub.tree,print.progress=F ) }, error=function(err) NA)
					#tryCatch( expr= {model <- two.b.pls(allo.Csize.bats[[anatomy[j]]][sub.tree$tip], pcoas[[ecology[k]]]$points[sub.tree$tip,],print.progress=F ) }, error=function(err) NA)

			} else  if( length(dim(allo.Csize.bats[[anatomy[j]]]) )<3 ) {
				tryCatch( expr= {model <- phylo.integration(allo.Csize.bats[[anatomy[j]]][sub.tree$tip,], pcoas[[ecology[k]]]$points[sub.tree$tip,], phy=sub.tree,print.progress=F ) }, error=function(err) NA)
					#tryCatch( expr= {model <- two.b.pls(allo.Csize.bats[[anatomy[j]]][sub.tree$tip,], pcoas[[ecology[k]]]$points[sub.tree$tip,],print.progress=F ) }, error=function(err) NA)				
			} else {
				tryCatch( expr= {model <- phylo.integration(allo.Csize.bats[[anatomy[j]]][,,sub.tree$tip], pcoas[[ecology[k]]]$points[sub.tree$tip,], phy=sub.tree,print.progress=F ) }, error=function(err) NA)
			}
			results.Z[[k]][i,j] <- model$Z[[1]] ; results.p[[k]][i,j] <- model$P.value[[1]]
		}
	}


	print(i)
	rm(allo.Csize)

}

library(reshape2) # 1.4.4
# Package for data curation
					 

melt.z.flight <- melt(results.Z[[2]])
melt.p.flight <- melt(results.p[[2]])
melt.z.flight <- cbind(melt.z.flight,rep('flight',dim(melt.z.flight)[1]))
colnames(melt.z.flight)[4]<-'eco'

melt.z.foot <- melt(results.Z[[3]])
melt.p.foot <- melt(results.p[[3]])
melt.z.foot <- cbind(melt.z.foot,rep('foot',dim(melt.z.foot)[1]))
colnames(melt.z.foot)[4]<-'eco'

# Produce curated data objects representing the strength and significance of bat anatomical adaptation to flight and roosting style variety. 

melt.z <- rbind(melt.z.flight,melt.z.foot)
melt.p <- rbind(melt.p.flight,melt.p.foot)

colour <- rep('black',length(melt.z$value))
colour[which(melt.p$value >0.05)] <- 'NA'
melt.z$value <- melt.z$value/10

colour<-colour[complete.cases(melt.z)]
melt.z<-melt.z[complete.cases(melt.z),]
melt.z <- cbind(melt.z,colour)


library(ggplot2) # 3.4.2
# Plotting package 
library(ggforce) # 0.4.1
# Package for constrained violin plots (if the user wishes to project individual data points into violin plots)

# The order in the violin plots we are about to produce needs to be 
# Gross wing proportions at top, then Handwing proportions, then handwing, radius, humerus, femur, tibia

library(ggnewscale) # 0.4.9
# Facilitates multiple overprinted axes in ggplots	 

library(jpeg)# version 0.1-10
# This package is for reading background graphics. These graphics are not provided in this GitHub repository, so the user seeking to replicate our results may dispense with this line. 

#img<- readJPEG('bat_bird_diagram2_crop.jpg')

# The user will not possess this image; they can dispense with running these lines

library(cowplot) # version 1.1.1
# Plot curation

library(magick) # version 2.8.2
# Graphics curation

blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(linewidth=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Blank background

n<-1
s<-1
p<-1

max.y<-1.6
sig <- quantile ( (n+(melt.z[-which(melt.z$eco=='flight'),]$value^p *s)[which( melt.z[-which(melt.z$eco=='flight'),]$colour=='black') ]), 0.05)

library(geomtextpath) # 0.1.1
# Facilitates curved text in ggplots

flight.significance <- rep('No', length(p.total$flight))
flight.significance[which(p.total$flight<0.05)]<- 'Yes'

roost.significance <- rep('No', length(p.total$roost))
roost.significance[which(p.total$roost<0.05)]<- 'Yes'
# Prepare significance data

melt.z<- cbind(melt.z, rep('NA', dim(melt.z)[1]))
colnames(melt.z)[6] <- 'significance'
melt.z$significance[which(melt.z$eco=='flight')] <- rep(flight.significance,each=100)
melt.z$significance[which(melt.z$eco=='foot')] <- rep(roost.significance,each=100)
# Prepare evolutionary covariance strength estimate data, ready to be presented as a violin plot representing ecological adaptation strength across the bat skeleton. 


bat.annot <-
ggplot(aes( x=Var2, y=((value^p) *s)+n , fill=eco,colour=Var2 ,linewidth=significance,linetype=significance), data=melt.z[which(melt.z$eco=='foot'),])+
geom_violin( trim=T  ,width=.75, position='dodge', linewidth=1, alpha=.5)+
scale_colour_manual( values = c("#0072B2","#0072B2","#0072B2","#CC79A7","#CC79A7","#0072B2","#0072B2"))+
#geom_text(aes(x=0.5,y=1,label='test'))+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.40,.40), y=c(sig-0.1,sig+0.1) )), aes(x=x,y=y),linetype='0',label='0.17')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.40,.40), y=c(1.39,1.41) )), aes(x=x,y=y),linetype='0',label='0.4')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.40,.40), y=c(.9,1.1) )), aes(x=x,y=y),linetype='0',label='0.0')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.15,.15), y=c(1.07,1.27) )), aes(x=x,y=y),linetype='0',label= 'paste(italic("Z/"), sqrt(n))', parse=T, size=6 )+

geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(1.45,1.45), y=c(sig-0.1,1.3) )), aes(x=x,y=y),linetype='0',label='Gross wing proportions')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(2.45,2.45), y=c(sig-0.1,1.4) )), aes(x=x,y=y),linetype='0',label='Handwing proportions')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(3.35,3.35), y=c(sig-0.4,1.4) )), aes(x=x,y=y),linetype='0',label='Handwing size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(4.35,4.35), y=c(sig-0.6,1.4) )), aes(x=x,y=y),linetype='0',label='Radius size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(5.40,5.40), y=c(sig-0.6,1.4) )), aes(x=x,y=y),linetype='0',label='Humerus size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(6.40,6.40), y=c(sig-0.3,1.4) )), aes(x=x,y=y),linetype='0',label='Femur size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(7.40,7.40), y=c(sig-0.3,1.4) )), aes(x=x,y=y),linetype='0',label='Tibia size')+
coord_polar(start = -1.57, direction=1, clip='on') +
#coord_curvedpolar(start = -1.57, direction=1, clip='on')+
scale_x_discrete(expand = c(1.2,1),limits=c('gross','Handwing.3D','handwing','radius','humerus','femur','tibia'), guide = guide_axis(angle = 90))+
scale_y_continuous(expand = c(0,0), limit=c(0,max.y) )+
scale_fill_manual(values = c('white','black'), labels = c("Flight-style", "Roosting/Foot-use"))+
geom_violin(aes( x=Var2, y=((value^p) *s)+n ,colour=Var2 ,linewidth=significance,linetype=significance), data=melt.z[which(melt.z$eco=='flight'),], trim=T  ,width=.75, position='dodge', linewidth=1, alpha=.75)+
scale_linetype_manual( values=c('dotted','solid'))+
annotate('segment', x = 0.5, xend = 7.5, y = sig, yend = sig, linewidth=3, colour = "black",alpha = 0.25, linetype='longdash')+
annotate('segment', x = 0.5, xend = 7.5, y = 1.4, yend = 1.4, linewidth=1, colour = "black",alpha = 0.25)+
annotate('segment', x = 0.5, xend = 7.5, y = 1, yend = 1, linewidth=1, colour = "black",alpha = 0.25)+
annotate('segment', x = 0.5, xend = 7.5, y = 1, yend = 1, linewidth=1, colour = "black",alpha = 0.25)+
guides(fill = guide_legend(override.aes = list(alpha = .25, shape=21)), colour = "none")+
labs(fill='Ecology', linetype='P < 0.05')+
  theme(
	legend.position='bottom',
	axis.text.x=element_blank(),
  	axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	plot.margin=unit(c(0,0,0,0), "pt"))
# Produce violin plot

took_legend <- get_legend(bat.annot)
# Extract legend


bat.annot <-
ggplot(aes( x=Var2, y=((value^p) *s)+n , fill=eco,colour=Var2 ,linewidth=significance,linetype=significance), data=melt.z[which(melt.z$eco=='foot'),])+
geom_violin( trim=T  ,width=.75, position='dodge', linewidth=1, alpha=.5)+
scale_colour_manual(values = c("#0072B2","#0072B2","#0072B2","#CC79A7","#CC79A7","#0072B2","#0072B2"))+
#geom_text(aes(x=0.5,y=1,label='test'))+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.40,.40), y=c(sig-0.1,sig+0.1) )), aes(x=x,y=y),linetype='0',label='0.17')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.40,.40), y=c(1.39,1.41) )), aes(x=x,y=y),linetype='0',label='0.4')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.40,.40), y=c(.9,1.1) )), aes(x=x,y=y),linetype='0',label='0.0')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(.15,.15), y=c(1.07,1.27) )), aes(x=x,y=y),linetype='0',label= 'paste(italic("Z/"), sqrt(n))', parse=T, size=6 )+

geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(1.45,1.45), y=c(sig-0.1,1.3) )), aes(x=x,y=y),linetype='0',label='Gross wing proportions')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(2.45,2.45), y=c(sig-0.1,1.4) )), aes(x=x,y=y),linetype='0',label='Handwing proportions')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(3.35,3.35), y=c(sig-0.4,1.4) )), aes(x=x,y=y),linetype='0',label='Handwing size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(4.35,4.35), y=c(sig-0.6,1.4) )), aes(x=x,y=y),linetype='0',label='Radius size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(5.40,5.40), y=c(sig-0.6,1.4) )), aes(x=x,y=y),linetype='0',label='Humerus size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(6.40,6.40), y=c(sig-0.3,1.4) )), aes(x=x,y=y),linetype='0',label='Femur size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(7.40,7.40), y=c(sig-0.3,1.4) )), aes(x=x,y=y),linetype='0',label='Tibia size')+
coord_polar(start = -1.57, direction=1, clip='on') +
#coord_curvedpolar(start = -1.57, direction=1, clip='on')+
scale_x_discrete(expand = c(1.2,1),limits=c('gross','Handwing.3D','handwing','radius','humerus','femur','tibia'), guide = guide_axis(angle = 90))+
scale_y_continuous(expand = c(0,0), limit=c(0,max.y) )+
scale_fill_manual(values = c('white','black'))+
geom_violin(aes( x=Var2, y=((value^p) *s)+n ,colour=Var2 ,linewidth=significance,linetype=significance), data=melt.z[which(melt.z$eco=='flight'),], trim=T  ,width=.75, position='dodge', linewidth=1, alpha=.75)+
scale_linetype_manual( values=c('dotted','solid'))+
annotate('segment', x = 0.5, xend = 7.5, y = sig, yend = sig, linewidth=3, colour = "black",alpha = 0.25, linetype='longdash')+
annotate('segment', x = 0.5, xend = 7.5, y = 1.4, yend = 1.4, linewidth=1, colour = "black",alpha = 0.25)+
annotate('segment', x = 0.5, xend = 7.5, y = 1, yend = 1, linewidth=1, colour = "black",alpha = 0.25)+
annotate('segment', x = 0.5, xend = 7.5, y = 1, yend = 1, linewidth=1, colour = "black",alpha = 0.25)+
  theme(
	legend.position='none',
	axis.text.x=element_blank(),
  	axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	plot.margin=unit(c(0,0,0,0), "pt"))
# plot sans legend





  ggdraw() +
  draw_plot(blank) +
#  draw_image(img,x=0.20,y=0.20,width=0.60,height=0.60)+
  draw_plot(bat.annot, x=-.10,y=-0.25,height=1.5, width=1.3)
# An example of the plot layout, with background graphics 


# The following lines will repeat the above in birds, and dense commenting will be desisted. 

# Load some requisite datasets
load('tree.22.10.2022.RData')
load('tree.names.22.10.2022.RData')
load('masses.22.10.2022.RData')
load('Csize.22.10.2022.RData')
load('coords.22.10.2022.RData')
# Done 

# Ensure names match the phylogeny. 
names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done

for(i in 1:length(GPA.coords)){
	dimnames(GPA.coords[[i]])[[3]]<- tree_names#dimnames(GPA.coords[[i]])
}


Csize.birds<-GPA.Csize[c('humerus','radius','carpometacarpus','femur','tibiotarsus')]

names(Csize.birds)[c(3,5)]<-c('handwing','tibia')



# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest. 


masses.birds <- masses
pruned.tree.birds <- pruned.tree

setwd()
bird.flight.styles <- read.csv('flight_masses_22_10_2022_plus_A.csv')
bird.foot.use <- read.csv('standardised_foot_scores_11_07_2023.csv')
bird.eco <- read.csv('Eco_meta_data.csv')


match.vector <- match(gsub('_.*','' ,tree_names),bird.eco$Prum_genus )
match.vector <- match.vector[complete.cases(match.vector)]

binary.flight.scores <- bird.flight.styles[,c(3:13)]
binary.flight.scores<-as.matrix(binary.flight.scores)
# Prepare the data as a matrix to make it easy to index. 
binary.flight.scores[which(binary.flight.scores=='?' | binary.flight.scores=='' )] <- 0 
binary.flight.scores[grep(binary.flight.scores,pattern='\\?')] <- 0 
binary.flight.scores <- apply(binary.flight.scores,2,FUN=as.numeric)
# I am making sure the data class is numeric
rownames(binary.flight.scores)<-bird.flight.styles$X
# Set the row names of the flight style matrix to our taxa 

binary.foot.scores <- bird.foot.use[,c(3:17)]
binary.foot.scores<-as.matrix(binary.foot.scores)
# Prepare the data as a matrix to make it easy to index. 
binary.foot.scores[which(binary.foot.scores=='?' | binary.foot.scores=='' | is.na(binary.foot.scores))] <- 0 
binary.foot.scores[grep(binary.foot.scores,pattern='\\?')] <- 0 

binary.foot.scores <- apply(binary.foot.scores,2,FUN=as.numeric)
# I am making sure the data class is numeric
rownames(binary.foot.scores)<-bird.foot.use$X


# The bird dataset has been loaded. 

anatomy <- c(names(Csize.birds),'carpometacarpus','brachial')

ecology <- c('flight','foot')



# First, we are going to compute an analysis with all taxa considered, in order guage the significance
# of relationships. 
# Thereafter, we will produce a confidence interval of Z values (effect size scores) by assessing 
# various subsamples. 

p.total <- list()
for(i in 1:length(ecology)){
	p.total[[i]] <- matrix(NA, 1, length(anatomy))
	colnames(p.total[[i]])<-anatomy
}
names(p.total)<-ecology


allo.Csize <- get.residual.Csize( array = Csize.birds, masses = masses, phylogeny = pruned.tree, taxa=names(masses)  )
allo.Csize.birds <- allo.Csize
brachial <-  Csize.birds$humerus[names(masses)]/GPA.Csize$ulna[names(masses)] # The 'brachial index' is a the humerus to ulna length ratio, and is widely perceived as an adaptive indicator of flight kinematics in birds 
carpometacarpus <-GPA.coords$carpometacarpus[,,names(masses)]
allo.Csize.birds[[6]] <- carpometacarpus
names(allo.Csize.birds)[6] <- c('carpometacarpus') 
allo.Csize.birds[[7]] <- brachial
names(allo.Csize.birds)[7] <- c('brachial') 



mat <- binary.flight.scores[names(masses),]
var<-apply(mat,2,sd)
mat <- mat[,which(var>0)]
flight.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
if(length(which(is.na(flight.dist)==T))>0){
	flight.dist[is.na(flight.dist)]<-0
}
flight.pcoa <- cmdscale(flight.dist, eig=T)
k<-max(which(flight.pcoa$eig/sum(flight.pcoa$eig)>=0.05))
flight.pcoa <- cmdscale( flight.dist, k=k  ,eig = T)


mat <- binary.foot.scores[names(masses),]
var<-apply(mat,2,sd)
mat <- mat[,which(var>0)]
foot.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
if(length(which(is.na(foot.dist)==T))>0){
	roost.dist[is.na(foot.dist)]<-0
}
foot.pcoa <- cmdscale(foot.dist, eig=T)
k<-max(which(foot.pcoa$eig/sum(foot.pcoa$eig)>=0.05))
foot.pcoa <- cmdscale( foot.dist, k=k  ,eig = T)

pcoas <- list(flight.pcoa,foot.pcoa)
names(pcoas)<-ecology

i<-1
for(j in 1: length(anatomy)){
	for(k in 1:length(ecology)){
		if( is.null(dim(allo.Csize.birds[[anatomy[j]]])) ){
			tryCatch( expr= {model <- phylo.integration(allo.Csize.birds[[anatomy[j]]][pruned.tree$tip], pcoas[[ecology[k]]]$points[pruned.tree$tip,], phy=pruned.tree,print.progress=F ) }, error=function(err) NA)
				#tryCatch( expr= {model <- two.b.pls(allo.Csize.birds[[anatomy[j]]][pruned.tree$tip], pcoas[[ecology[k]]]$points[pruned.tree$tip,],print.progress=F ) }, error=function(err) NA)

		} else  if( length(dim(allo.Csize.birds[[anatomy[j]]]) )<3 ) {
			tryCatch( expr= {model <- phylo.integration(allo.Csize.birds[[anatomy[j]]][pruned.tree.birds$tip,], pcoas[[ecology[k]]]$points[pruned.tree$tip,], phy=pruned.tree,print.progress=F ) }, error=function(err) NA)
				#tryCatch( expr= {model <- two.b.pls(allo.Csize.birds[[anatomy[j]]][pruned.tree.birds$tip,], pcoas[[ecology[k]]]$points[pruned.tree$tip,],print.progress=F ) }, error=function(err) NA)				
		} else {
			tryCatch( expr= {model <- phylo.integration(allo.Csize.birds[[anatomy[j]]][,,pruned.tree$tip], pcoas[[ecology[k]]]$points[pruned.tree$tip,], phy=pruned.tree,print.progress=F ) }, error=function(err) NA)
		}
		 p.total[[k]][i,j] <- model$P.value[[1]] #results.Z[[k]][i,j] <- model$Z[[1]] ;
	}
}

# p.total demonstrates that the gross proportions of the bat wing are consistently significantly associated to all ecological categories. 
# The carpometacarpus size is relevant to flight style, as is the femur and the carpometacarpus shape 
# The femur size and tibia size are relevant to foot-use ecology 





results.Z <- list()
for(i in 1:length(ecology)){
	results.Z[[i]] <- matrix(NA, 100, length(anatomy))
	colnames(results.Z[[i]])<-anatomy
}
names(results.Z)<-ecology

results.p <- results.Z 


for(i in 1: 100){


	a<-1
	while(a<2){
		subsample <- sample(tree_names,100)
		
		tryCatch( expr= {allo.Csize <- get.residual.Csize( array = Csize.birds, masses = masses, phylogeny = pruned.tree, taxa=subsample  ) }, error=function(err) NA)
		if(exists('allo.Csize')){
			a<-2
		}
		sub.tree <- keep.tip(pruned.tree, subsample)
		
	}


	brachial <- Csize.birds$humerus[subsample]/GPA.Csize$ulna[subsample]


	carpometacarpus <-  GPA.coords$carpometacarpus[,,subsample]


	allo.Csize.birds <- allo.Csize	
	allo.Csize.birds[[6]] <- carpometacarpus
	allo.Csize.birds[[7]] <- brachial 
	names(allo.Csize.birds)[6:7] <- c('carpometacarpus','brachial') 
			 
	mat <- binary.flight.scores[subsample,]
	var<-apply(mat,2,sd)
	mat <- mat[,which(var>0)]
	flight.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
	if(length(which(is.na(flight.dist)==T))>0){
		flight.dist[is.na(flight.dist)]<-0
	}
	flight.pcoa <- cmdscale(flight.dist, eig=T)
	k<-max(which(flight.pcoa$eig/sum(flight.pcoa$eig)>=0.05))
	flight.pcoa <- cmdscale( flight.dist, k=k  ,eig = T)


	mat <- binary.foot.scores[subsample,]
	var<-apply(mat,2,sd)
	mat <- mat[,which(var>0)]
	foot.dist <- daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
	if(length(which(is.na(foot.dist)==T))>0){
		foot.dist[is.na(foot.dist)]<-0
	}
	foot.pcoa <- cmdscale(foot.dist, eig=T)
	k<-max(which(foot.pcoa$eig/sum(foot.pcoa$eig)>=0.05))
	foot.pcoa <- cmdscale( foot.dist, k=k  ,eig = T)

	pcoas <- list(flight.pcoa,foot.pcoa)
	names(pcoas)<-ecology

	for(j in 1: length(anatomy)){
		for(k in 1:length(ecology)){
			if( is.null(dim(allo.Csize.birds[[anatomy[j]]])) ){
				tryCatch( expr= {model <- phylo.integration(allo.Csize.birds[[anatomy[j]]][sub.tree$tip], pcoas[[ecology[k]]]$points[sub.tree$tip,], phy=sub.tree,print.progress=F ) }, error=function(err) NA)
					#tryCatch( expr= {model <- two.b.pls(allo.Csize.birds[[anatomy[j]]][sub.tree$tip], pcoas[[ecology[k]]]$points[sub.tree$tip,],print.progress=F ) }, error=function(err) NA)

			} else if( length(dim(allo.Csize.birds[[anatomy[j]]]) )<3 ) {
				tryCatch( expr= {model <- phylo.integration(allo.Csize.birds[[anatomy[j]]][sub.tree$tip,], pcoas[[ecology[k]]]$points[sub.tree$tip,], phy=sub.tree,print.progress=F ) }, error=function(err) NA)
					#tryCatch( expr= {model <- two.b.pls(allo.Csize.birds[[anatomy[j]]][sub.tree$tip,], pcoas[[ecology[k]]]$points[sub.tree$tip,],print.progress=F ) }, error=function(err) NA)				
			} else {
				tryCatch( expr= {model <- phylo.integration(allo.Csize.birds[[anatomy[j]]][,,sub.tree$tip], pcoas[[ecology[k]]]$points[sub.tree$tip,], phy=sub.tree,print.progress=F ) }, error=function(err) NA)

			}
			results.Z[[k]][i,j] <- model$Z[[1]] ; results.p[[k]][i,j] <- model$P.value[[1]]
		}
	}


	print(i)
	rm(allo.Csize)
}




melt.z.flight <- melt(results.Z[[2]])
melt.p.flight <- melt(results.p[[2]])
melt.z.flight <- cbind(melt.z.flight,rep('flight',dim(melt.z.flight)[1]))
colnames(melt.z.flight)[4]<-'eco'

melt.z.foot <- melt(results.Z[[3]])
melt.p.foot <- melt(results.p[[3]])
melt.z.foot <- cbind(melt.z.foot,rep('foot',dim(melt.z.foot)[1]))
colnames(melt.z.foot)[4]<-'eco'

melt.z <- rbind(melt.z.flight,melt.z.foot)
melt.p <- rbind(melt.p.flight,melt.p.foot)

colour <- rep('black',length(melt.z$value))
colour[which(melt.p$value >0.05)] <- 'NA'
melt.z$value <- melt.z$value/10
#melt.z$value[which(melt.z$value<0)]<-NA

colour<-colour[complete.cases(melt.z)]
melt.z<-melt.z[complete.cases(melt.z),]


flight.significance <- rep('No', length(p.total$flight))
flight.significance[which(p.total$flight<0.05)]<- 'Yes'

foot.significance <- rep('No', length(p.total$foot))
foot.significance[which(p.total$foot<0.05)]<- 'Yes'

melt.z<- cbind(melt.z, rep('NA', dim(melt.z)[1]))
colnames(melt.z)[5] <- 'significance'
melt.z$significance[which(melt.z$eco=='flight')] <- rep(flight.significance,each=100)
melt.z$significance[which(melt.z$eco=='foot')] <- rep(foot.significance,each=100)

bird.annot <-
ggplot(aes( x=Var2, y=((value^p) *s)+n , fill=eco,colour=Var2,linewidth=significance,linetype=significance ), data=melt.z[which(melt.z$eco=='foot'),])+
geom_violin( trim=T  ,width=.75, position='dodge', linewidth=1, alpha=.5)+
scale_colour_manual(values = c("#0072B2","#0072B2","#0072B2","#CC79A7","#CC79A7","#0072B2","#0072B2"))+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(7.55,7.55), y=c(sig-0.1,sig+0.1) )), aes(x=x,y=y),linetype='0',label='0.17')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(7.55,7.55), y=c(1.39,1.41) )), aes(x=x,y=y),linetype='0',label='0.4')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(7.55,7.55), y=c(.9,1.1) )), aes(x=x,y=y),linetype='0',label='0.0')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(7.75,7.75), y=c(1.07,1.27) )), aes(x=x,y=y),linetype='0',label= 'paste(italic("Z/"), sqrt(n))', parse=T , size=6)+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(0.55,0.55), y=c(sig-0.1,1.3) )), aes(x=x,y=y),linetype='0',label='Tibia size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(1.55,1.55), y=c(sig-0.1,1.4) )), aes(x=x,y=y),linetype='0',label='Femur size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(2.55,2.55), y=c(sig-0.5,1.4) )), aes(x=x,y=y),linetype='0',label='Humerus size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(3.65,3.65), y=c(sig-0.3,1.4) )), aes(x=x,y=y),linetype='0',label='Radius size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(4.60,4.60), y=c(sig-0.2,1.4) )), aes(x=x,y=y),linetype='0',label='Carpometacarpus size')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(5.55,5.55), y=c(sig-0.1,1.4) )), aes(x=x,y=y),linetype='0',label='Carpometacarpus shape')+
geom_textpath(inherit.aes=F, data= data.frame(cbind(x=c(6.55,6.55), y=c(sig-0.1,1.4) )), aes(x=x,y=y),linetype='0',label='Brachial index')+
coord_polar(start = 1.57, direction=1, clip='on') +
scale_x_discrete(expand = c(1.2,1),limits=rev(c('brachial','carpometacarpus','handwing','radius','humerus','femur','tibia')), guide = guide_axis(angle = 90))+
scale_y_continuous(expand = c(0,0), limit=c(0,max.y) )+
scale_fill_manual(values = c('white','black'))+
geom_violin(aes( x=Var2, y=((value^p) *s)+n ,colour=Var2,linewidth=significance,linetype=significance), data=melt.z[which(melt.z$eco=='flight'),], trim=T  ,width=.75, position='dodge', linewidth=1, alpha=.75)+
scale_linetype_manual( values=c('dotted','solid'))+
annotate('segment', x = 0.5, xend = 7.5, y = sig, yend = sig, linewidth=3, colour = "black",alpha = 0.25, linetype='longdash')+
annotate('segment', x = 0.5, xend = 7.5, y = 1.4, yend = 1.4, linewidth=1, colour = "black",alpha = 0.25)+
annotate('segment', x = 0.5, xend = 7.5, y = 1, yend = 1, linewidth=1, colour = "black",alpha = 0.25)+
annotate('segment', x = 0.5, xend = 7.5, y = 1, yend = 1, linewidth=1, colour = "black",alpha = 0.25)+
  theme(
	legend.position='none',
	axis.text.x=element_blank(),
  	axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	plot.margin=unit(c(0,0,0,0), "pt"))
# Done




dev.new(width=18,height=13,unit='cm')
# Produce image

  ggdraw() +
  draw_plot(blank) +
	draw_plot(took_legend,width=1,height=1.9)+
 # draw_image(img,x=0.2,y=-0.00,width=0.58,height=.9)+
  draw_plot(bat.annot, x=0.03,y=-0.49,height=2,width=1.1)+
  draw_plot(bird.annot, x=-0.13,y=-0.49,height=2,width=1.1)+
geom_text(data=data.frame(cbind(x=c(0.1,0.9),y=c(0.95,0.95))), aes(x=x,y=y), label=c('a','b'), size = 30/.pt, color = "black", fontface = "bold")
# Combine bat and ecological adaptation plots and legend onto a single figure background

ggsave(file='polar_plot_04_05_2024.pdf',height=26,width=36,unit='cm',dpi=300)
# Save the figure 


# Script concludes 


