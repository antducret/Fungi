
# PREPARE SESSION ---------------------------------------------------------

# Set directory ?
# Librairies ?
setwd(dir = "/home/filsdufrere/Documents/Multivariate Statistics in R/datasets_fungi/Fungi/Projet_Rstudio")
library(tidyverse)

world_map = map_data("world")


# PATH ---------------------------------------------------------------
PATH_SPE = "../Data/occurence_tot.csv"
PATH_SPE_BENTHIC = "../Data/occurence_benthic.csv"
PATH_SPE_PELAGIC = "../Data/occurence_pelagic.csv"
PATH_ENV = "../Data/env.csv"

# Load data ---------------------------------------------------------------
  # Environmental parameters
  env <- read.csv(PATH_ENV, sep = ',', skip = 1 )
  env <- env[2205:2214]   # Keep only environmental parameters
  env 
  
  #speurrences
  
    # Total
    spe <- read.csv(PATH_SPE, sep = ',', row.names = NULL)
    name_species <- spe[,2:11]
    
    spe <- spe[spe$Class != '',]
    spe$Class <- as.factor(spe$Class)
    spe <- aggregate(spe[,-(0:11)], by = list(spe$Class), FUN = "sum")
    site_names <- colnames(spe[,-1], prefix = "")
    namecol <- append(1:56, "Tot", 0)
    colnames(spe) <- append(namecol, "Class", 0) # Pourquoi ?
    
    head(spe)
    

#asp = 1, 
#type = "p", 
#main = "Site Locations", 
#xlab = "° longitude", 
#
#ylim = c(-90,90),
#xlim = c(-180,180),
# Data exploration SPE  ---------------------------------------------------------------
    
    # Summary
      summary(spe)
      str(spe)
  
    # Proportion of zeros in the community data set
      sum(spe == 0) / (nrow(spe) * ncol(spe))
      
    # relative presence of each species
            
        # total presence for each species
          spe.pres <- apply(spe[, 2:58]> 0, 1, sum)
          sites.pres <- apply(spe[, 2:58]> 0, 2, sum)
        # Sort the results in increasing order
          sort(spe.pres)
        
        # Plot histograms
        barplot(sites.pres[-(1)], 
               las = 1,
               xlab = "Sites",
               ylab = "Species richness",
               col = gray(5 : 0 / 5),
               horiz=F,
       )
        # Data exploration ENV  ---------------------------------------------------------------
        # Summary 
        summary(env)
        str(env) # asfactor pour N/Y sediment ?  
        
        # Map sites (TODO : Improve)
        coo = env[,(1:2)]
        map <- ggplot(world_map, aes(x = long, y = lat, group = group))  +
          geom_path()+
          scale_y_continuous(breaks = (-2:2) * 30) +
          scale_x_continuous(breaks = (-4:4) * 45) +
          geom_point(data = coo, aes(x=Longitude..degrees., y=Latitide..degrees., size = sites.pres[-1],colour = "red"), inherit.aes = FALSE) +
          xlab("° Longitude") + ylab("° Latitude") +
          labs(colour ='Sites of study',size = "Number of occurrences of classes" , title = "Location of sites")
        map 
        
         

#Supervised classification
    #Classification trees, Random Forest 
    # Hierarchical (dendrograms) or not (k-means)

spe.norm = decostand(spe[,2:58], 'normalize')

spe.hel  = vegdist(spe.norm, 'hel') ## psk bcp 0

# Compute and plot complete-linkage agglomerative clustering

spe.hel.complete <- hclust(spe.hel, method = "complete")
plot(spe.hel.complete, main = "Chord - Complete linkage")
spe.hel.complete

# Compute UPGMA clustering
spe.hel.UPGMA <- hclust(spe.hel, method = "average")
plot(spe.hel.UPGMA, main = "Chord - UPGMA")


# Compute UPGMC clustering
spe.hel.centroid <- hclust(spe.hel, method = "centroid")
plot(spe.hel.centroid, main = "Chord - Centroid")


# Compute Ward's minimum variance clustering
spe.hel.ward <- hclust(spe.hel, method = "ward.D2")
plot(spe.hel.ward,  main = "Chord - Ward")


#k-means clustering
spe.hel.k<-  kmeans(spe.hel, centers=5)
spe.hel.k


# Cophenetic correlations -------------------------------------------------
# Single linkage clustering
spe.hel.single.coph <- cophenetic(spe.hel.single)
cor(spe.hel, spe.hel.single.coph)


# Complete linkage clustering
spe.hel.comp.coph <- cophenetic(spe.hel.complete)
cor(spe.hel, spe.hel.comp.coph)


# Average clustering
spe.hel.UPGMA.coph <- cophenetic(spe.hel.UPGMA)
cor(spe.hel, spe.hel.UPGMA.coph)


# Ward clustering
spe.hel.ward.coph <- cophenetic(spe.hel.ward)
cor(spe.hel, spe.hel.ward.coph)

library(NbClust)


Nb.complete <-NbClust(spe[,2:58], diss=spe.hel, distance = NULL, min.nc=1, max.nc=16, 
                   method = "average", index = "ch")

Nb.complete
plot(Nb.complete$All.index, xlab="number of clusters", ylab="Calinski and Harabasz index")

single.dend <- as.dendrogram(spe.hel.complete)
plot(single.dend)
#Unconstrained ordination
    #CA
      
#Constrained ordination 
    # CCA
        
            
