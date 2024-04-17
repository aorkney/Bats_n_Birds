# The purpose of this script is to 
# produce a phenogram ancestral state reconstruction plot of bird and bat dispersion over
# wing:(wing+leg), and colour the branches by the inferred rate changes


library( geomorph ) # 4.0.5
library( ape ) # 5.7-1
library( ggplot2 ) # 3.4.1
library( ggdendro ) # 0.1.23
library( dendextend ) # 1.17.1
library(zoo) # 1.8-12
library(dplyr) # 1.1.1
library( phytools ) # 1.5-1
library( mvMORPH) # 1.1.8
library(cowplot) # 1.1.1
library(geiger) # 2.0.10
library(caTools) # 1.18.2
library(ggplot2) # 3.4.1

# Load required packages

# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function


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
# Coordinates and centroid sizes have been extrated

setwd()
# Set work directory to the location of the phylogeny by Shi and Rabosky 

bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# This phylogeny far exceeds the number of taxa for which we have landmark constellations
# we must therefore prune the phylogeny  

pruned.tree.bats <- keep.tip(bat.tree,taxa)






bat.proportions <- cbind((GPA.Csize$humerus)+(GPA.Csize$radius)+(GPA.Csize$handwing),
(GPA.Csize$femur)+(GPA.Csize$tibia))
bat.proportions <- bat.proportions/rowSums(bat.proportions)



# Set the work directory
setwd()
# The user will need to customise this


# Load some requisite datasets
load('tree.22.10.2022.RData')
load('tree.names.22.10.2022.RData')
load('Csize.22.10.2022.RData')
# Done 

# Ensure names match the phylogeny. 
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done

Csize.birds<-GPA.Csize

avian.proportions <- cbind((Csize.birds$humerus)+(Csize.birds$radius)+(Csize.birds$carpometacarpus),
(Csize.birds$femur)+(Csize.birds$tibiotarsus))
avian.proportions <- avian.proportions/rowSums(avian.proportions)





make.dendr <- function(tree, trait, model){
	dendr <- dendro_data(as.dendrogram(force.ultrametric(tree)))
	lab.dat <- dendr$labels
	if(model=='BM'){
		fit <- fastAnc(tree,trait)
		Anc<- rep(NA, length(dendr$segments[,1]))
		Anc [which(dendr$segments$x==dendr$segments$xend)] <- fit[match(tree$edge[,1],names(fit))]
		Anc [which(dendr$segments$x==dendr$segments$xend)-1] <- fit[match(tree$edge[,1],names(fit))]

		dendr.mod<-dendr$segments/2
		dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)-1] <- fit[match(tree$edge[,1],names(fit))] # This defines the nodal positions of the nodes
		dendr.mod$xend[which(dendr$segments$x==dendr$segments$xend)-1] <- dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)-1]# This collapses nodes by defining their ends as the same as their starts
		dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)] <-   fit[match(tree$edge[,1],names(fit))] # This defines the nodal positions of the start of each branch
		dendr.mod$xend[which(dendr$segments$yend==0)] <- trait[match(lab.dat$label,names(trait)) ] # This defines the final position of each terminal branch
	}

	if(model=='EB'){
		mod <- mvEB(tree,trait)
		fit<-estim(tree, data= trait, object=mod, asr=TRUE)

		Anc <- rep(NA, length(dendr$segments[,1])) 
		Anc [which(dendr$segments$x==dendr$segments$xend)] <- fit$estim[match(tree$edge[,1],rownames(fit$estim))] 
		Anc [which(dendr$segments$x==dendr$segments$xend)-1] <- fit$estim[match(tree$edge[,1],rownames(fit$estim))]

		dendr.mod<-dendr$segments/2
		dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)-1] <- fit$estim[match(tree$edge[,1],rownames(fit$estim))] # This defines the nodal positions of the nodes
		dendr.mod$xend[which(dendr$segments$x==dendr$segments$xend)-1] <- dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)-1]# This collapses nodes by defining their ends as the same as their starts
		dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)] <-   fit$estim[match(tree$edge[,1],rownames(fit$estim))] # This defines the nodal positions of the start of each branch
		dendr.mod$xend[which(dendr$segments$yend==0)] <- trait[match(lab.dat$label,names(trait)) ] # This defines the final position of each terminal branch
	
	}


	for(i in 1:dim(dendr.mod)[1]){
		if(dendr$segments[i,]$yend !=0 ){ # We do not want to move distal tips of branches that hit zero
			if( is.na(match(dendr.mod[i,]$xend,dendr.mod$x))==T  ){ # We check that the branch has no descendants 
				choice <- which(dendr$segments$x==dendr$segments[i,]$xend) # We find a potential site it should attach to 
				target <- choice[which(dendr$segments[choice,]$y-dendr$segments[choice,]$yend == 0)][1] # We have found a node to target
				dendr.mod$xend[i]<- dendr.mod$x[target]
			
			}
		}	
	}

	Rate <- abs((dendr.mod$xend-dendr.mod$x))/(dendr.mod$y-dendr.mod$yend)

	return( cbind(dendr.mod, Anc,Rate ) )

}

dendr.bat <- make.dendr(tree=pruned.tree.bats, trait=bat.proportions[pruned.tree.bats$tip,1], model='EB')
dendr.bird <- make.dendr(tree=pruned.tree, trait=avian.proportions[pruned.tree$tip,1], model='BM')


dendr.col <- rbind(dendr.bat,dendr.bird)


plot.with.legend<-
ggplot()+
geom_segment(data=dendr.col[order(dendr.col$Rate),], aes(x = -y, y = x, xend = -yend, yend =xend ),linewidth=3, col='black' )+
geom_segment(data=dendr.col[order(dendr.col$Rate),], aes(x = -y, y = x, xend = -yend, yend =xend, col=Rate^.5 ),linewidth=1.5 )+
scale_color_gradientn(colors=rev(c('white','yellow','red','blue','darkblue')))+
#labs(col= expression(italic(sqrt(Rate))))+
labs(col= expression(italic(sigma)), x=expression(italic('Mya')), y=expression(italic('wing/(leg+wing)')))+
	theme(legend.position='right',
         panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.title= element_text(size=40,color='black'),
		legend.text= element_text(size=16,color='black'),
		legend.key.width=unit( 0.8,'cm'),
		legend.key.height=unit( 1,'cm'),
		axis.title.x=element_text(size=30, color='black'),
		axis.title.y=element_text(size=30, color='black'),

		#axis.title.y=element_blank(),
		#axis.text.x=element_blank(),
		#axis.text.y=element_blank(),
		#axis.ticks.y=element_blank(),
		#axis.ticks.x=element_blank(),
       )
plot.with.legend<-
plot.with.legend+
 guides(col = guide_colourbar(ticks.colour= 'black', ticks.linewidth=6/.pt, frame.colour = 'black',
  frame.linewidth = 4/.pt))

legend<- cowplot::get_legend(plot.with.legend)


plot.without.legend<-
ggplot()+
geom_segment(data=dendr.col[order(dendr.col$Rate),], aes(x = -y, y = x, xend = -yend, yend =xend ),linewidth=3.5, col='black' )+
geom_segment(data=dendr.col[order(dendr.col$Rate),], aes(x = -y, y = x, xend = -yend, yend =xend, col=Rate^.5  ),linewidth=1.5 )+
scale_color_gradientn(colors=rev(c('white','yellow','red','blue','darkblue')))+
labs(col= expression(italic(sigma)), x=expression(italic('Mya')), y=expression(italic('wing/(leg+wing)')))+
geom_text(aes(x=c(-40,-40),y=c(0.4,0.85),label=c('birds','bats')),color='black',size=10 ) +
	theme(legend.position='none',
         panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.title= element_text(size=40,color='white'),
		legend.text= element_text(size=16,color='white'),
		legend.key.width=unit( 0.8,'cm'),
		legend.key.height=unit( 1,'cm'),
		axis.line.x=element_line(color='black',linewidth=3/2),
		axis.line.y=element_line(color='black',linewidth=3/2),
		axis.ticks.x=element_line(color='black',linewidth=2),
		axis.ticks.y=element_line(color='black',linewidth=2),
		axis.title.x=element_text(size=30, color='black'),
		axis.title.y=element_text(size=30, color='black'),
		axis.text.x=element_text(size=16,color='black'),
		axis.text.y=element_text(size=16,color='black'),
		#axis.title.y=element_blank(),
		#axis.text.x=element_blank(),
		#axis.text.y=element_blank(),
		#axis.ticks.y=element_blank(),
		#axis.ticks.x=element_blank(),
       )

tree.plot<-
plot.without.legend+
draw_plot(legend , x = -72, y = 0.6, width = 15, height = 0.2 )

resample.dtt <- function(phy, data, sample.size, repeats){
	return <- list()
	
	times <- list()
	real <- list()
	sim <- list()
	for(i in 1:repeats){
		sample <- sample(phy$tip,sample.size)
		DTT <- dtt(phy=keep.tip(phy,sample) , data= data[sample] , nsim=1, plot=F)
		times[[i]] <- DTT$times
		real[[i]] <- DTT$dtt
		sim[[i]] <- DTT$sim	
	}
	
	return[[1]] <- unlist(times)
	return[[2]] <- unlist(real)
	return[[3]] <- unlist(sim)
	names(return)<-c('times','real','sim')
	return(return)
} 



DTT.bat <- resample.dtt(phy=pruned.tree.bats,data=bat.proportions[,1],sample.size=100,repeats=100)
DTT.bird <- resample.dtt(phy=pruned.tree,data=avian.proportions[,1],sample.size=100,repeats=100)

quant.bat <- runquantile(x=DTT.bat$sim[rev(order(DTT.bat$times))], probs=c(0.023,0.159,0.841,0.977), k=20)
quant.bird <- runquantile(x=DTT.bird$sim[rev(order(DTT.bird$times))], probs=c(0.023,0.159,0.841,0.977), k=20)

dtt.plot.data <- function(DTT.data, probs, k, span, sim){
	if(sim==T){
		quant <- runquantile(x=DTT.data$sim[rev(order(DTT.data$times))],probs=probs,k=k)
	} else {
		quant <- runquantile(x=DTT.data$real[rev(order(DTT.data$times))],probs=probs,k=k)
	}
	for(i in 1:length(probs)){
		if(i==1){
			fit <- loess(quant[,i] ~ DTT.data$times[rev(order(DTT.data$times))],span=span)
			smoothed <- cbind(fit$x,fit$fitted)
		} else {
			smoothed <- cbind(smoothed, loess(quant[,i] ~ DTT.data$times[rev(order(DTT.data$times))],span=span)$fitted)
		}
	}
	colnames(smoothed) <- c('times',probs)
	return(data.frame(smoothed))
}


ribbon.bat.sim <- dtt.plot.data(DTT.data=DTT.bat,probs=c(0.023,0.159,0.841,0.977), k=20,span=0.15, sim=T)

ribbon.bird.sim <- dtt.plot.data(DTT.data=DTT.bird,probs=c(0.023,0.159,0.841,0.977), k=20,span=0.15, sim=T)

ribbon.bat.real <- dtt.plot.data(DTT.data=DTT.bat,probs=c(0.023,0.159,0.841,0.977), k=20,span=0.15, sim=F)

ribbon.bird.real <- dtt.plot.data(DTT.data=DTT.bird,probs=c(0.023,0.159,0.841,0.977), k=20,span=0.15, sim=F)


library(png); library(grid)

setwd()

bird.img = readPNG('Bubo blakistoni_512x345.png')
bird.img[,,2][which(bird.img[,,2]==0)]<-NA
bird.img[,,2][which(bird.img[,,2]==1)]<-0
bird.img[,,2][which(bird.img[,,2]!=0)]<-NA


bat.img = readPNG('Vespertilionidae_Corynorhinus townsendii_test.png')
bat.img[,,2][which(bat.img[,,2]==0)]<-NA
bat.img[,,2][which(bat.img[,,2]==1)]<-0
bat.img[,,2][which(bat.img[,,2]!=0)]<-NA

#bat.g =  rasterGrob(bat.img[,,2], interpolate=F)
#bird.g =  rasterGrob(bird.img[,,2], interpolate=F)

bat.raster <- as.raster(bat.img[,,2], interpolate=F)
bird.raster <- as.raster(bird.img[,,2], interpolate=F)

back.bat <- cbind(ribbon.bat.sim, rep(mean(dendr.bat$Rate^.5,na.rm=T),dim(ribbon.bat.sim)[1] ))
colnames(back.bat)[6]<-'f'

bat.disp<-
ggplot()+
geom_ribbon(data= back.bat, aes(x=times, ymin= X0.023, ymax=X0.977, fill=f ) )+
scale_fill_gradientn( colors=rev(c('white','yellow','red','blue','darkblue')), limits=range(dendr.col$Rate[order(dendr.col$Rate)]^.5,na.rm=T) )+
geom_ribbon(data= ribbon.bat.real, aes(x=times, ymin= X0.023, ymax=X0.977), fill='black',alpha=.75)+
labs(x='',y='')+
lims(y=c(-0.1,1.8))+
scale_x_continuous(breaks=c(0,.25,.5,.75,1), limits=c(0,1))+
	theme(legend.position='none',
         panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.title= element_text(size=40,color='black'),
		legend.text= element_text(size=16,color='black'),
		legend.key.width=unit( 0.8,'cm'),
		legend.key.height=unit( 1,'cm'),
		axis.line.x=element_line(color='black',linewidth=3/2),
		axis.line.y=element_line(color='black',linewidth=3/2),
		axis.ticks.x=element_line(color='black',linewidth=2),
		axis.ticks.y=element_line(color='black',linewidth=2),
		axis.title.x=element_text(size=30, color='black'),
		axis.title.y=element_text(size=30, color='black'),
		axis.text.x=element_text(size=16,color='black'),
		axis.text.y=element_text(size=16,color='black')
       )


g <- ggplot_build(bat.disp)
bat.raster[bat.raster == "#000000"] <- unique(g$data[[1]]["fill"])[[1]][1]
bat.g =  rasterGrob(bat.raster, interpolate=F)

bat.disp <- bat.disp +annotation_custom(grob=bat.g, xmin=.45, xmax=.85, ymin=0.9, ymax=1.5)


back.bird <- cbind(ribbon.bird.sim, rep(mean(dendr.bird$Rate^.5,na.rm=T),dim(ribbon.bird.sim)[1] ))
colnames(back.bird)[6]<-'f'

bird.disp<-
ggplot()+
geom_ribbon(data= back.bird, aes(x=times, ymin= X0.023, ymax=X0.977, fill=f ) )+
scale_fill_gradientn( colors=rev(c('white','yellow','red','blue','darkblue')), limits=range(dendr.col$Rate[order(dendr.col$Rate)]^.5,na.rm=T) )+
geom_ribbon(data= ribbon.bird.real, aes(x=times, ymin= X0.023, ymax=X0.977), fill='black',alpha=.75)+
labs(x='',y='')+
lims(y=c(-0.1,1.8))+
scale_x_continuous(breaks=c(0,.25,.5,.75,1), limits=c(0,1))+
	theme(legend.position='none',
         panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.title= element_text(size=40,color='white'),
		legend.text= element_text(size=16,color='white'),
		legend.key.width=unit( 0.8,'cm'),
		legend.key.height=unit( 1,'cm'),
		axis.line.x=element_line(color='black',linewidth=3/2),
		axis.line.y=element_line(color='black',linewidth=3/2),
		axis.ticks.x=element_line(color='black',linewidth=2),
		axis.ticks.y=element_line(color='black',linewidth=2),
		axis.title.x=element_text(size=30, color='black'),
		axis.title.y=element_text(size=30, color='black'),
		axis.text.x=element_text(size=16,color='black'),
		axis.text.y=element_text(size=16,color='black')
       )


g <- ggplot_build(bird.disp)
bird.raster[bird.raster == "#000000"] <- unique(g$data[[1]]["fill"])[[1]][1]
bird.g =  rasterGrob(bird.raster, interpolate=F)

bird.disp <- bird.disp +annotation_custom(grob=bird.g, xmin=.3, xmax=.75, ymin=0.4, ymax=1.8)


library(ggpubr)

disp.plot <- ggarrange(bat.disp,bird.disp, ncol=1)
disp.plot <- annotate_figure(disp.plot, left = text_grob(expression(italic("subclade disparity")), 
              color = "black", face = "bold", size = 30, rot=90, vjust=1.75))
disp.plot <- annotate_figure(disp.plot, bottom = text_grob(expression(italic("relative age")), 
               color = "black", face = "bold", size = 30, vjust=-.75, hjust=0.3))

assembled <- ggarrange(tree.plot,disp.plot, labels=c('a','b'),font.label = list(size = 30, color = "black", face = "bold", family = NULL) )
#assembled <- assembled +bgcolor('white')


dev.new(height=13.5,width=20,unit='cm') 
assembled

>> Sub-sample to the square root of the overall number of species 


setwd()
ggsave(filename='Figure_pheno_02_20_2024_pale.pdf', width = 40,
  height = 27,
  units = c( "cm"),
  dpi = 300)
