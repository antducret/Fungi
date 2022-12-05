
# PREPARE SESSION ---------------------------------------------------------

# Set directory ?
# Librairies ?
setwd(dir = "/home/filsdufrere/Documents/Multivariate Statistics in R/datasets_fungi/Fungi/Projet_Rstudio")
#setwd("/home/anthoney/Documents/Master/MiR/Fungi/Projet_Rstudio/")
library(tidyverse)
library(vegan)

world_map = map_data("world")


# PATH ---------------------------------------------------------------
PATH_SPE = "../Data/occurence_tot.csv"
PATH_SPE_BENTHIC = "../Data/occurence_benthic.csv"
PATH_SPE_PELAGIC = "../Data/occurence_pelagic.csv"
PATH_ENV = "../Data/env.csv"

# Load & Clean data ---------------------------------------------------------------
  # ENV
    env <- read.csv(PATH_ENV, sep = ',', skip = 1 )
    env <- env[2205:2214]   # Keep only environmental parameters
    env$sediment <- as.factor(env$sediment)
    env 
  
  # SPE
    spe <- read.csv(PATH_SPE, sep = ',', row.names = NULL)
    name_species <- spe[,2:11]
    
    spe <- spe[spe$Class != '',]
    spe$Class <- as.factor(spe$Class)
    
    spe <- aggregate(spe[,-(0:11)], by = list(spe$Class), FUN = "sum")
    site_names <- colnames(spe[,-1], prefix = "")
    colnames(spe) <- append(1:56, c("Class","Tot"), 0)

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
  ## Hierarchical agglomerative clustering of the species abundance 
        
        # Compute Helinger distance
        spe.norm <- decostand(spe[,2:58], "normalize")
        spe.hel <- vegdist(spe.norm, "hel")
        
        # METHOD
          # CREATE DENDROGRAM
          # PLOT DENDROGRAM
          # COMPUTE COPHENETIC MATRIX
          # COMPUTE CORRELATION
        
        # Single method 
        spe.hel.single <- hclust(spe.hel, method = "single")
        plot(spe.hel.single, main = "Helinger - Single linkage")
        spe.hel.single.coph <- cophenetic(spe.hel.single)
        cor(spe.hel, spe.hel.single.coph)
        
        # Complete method 
        spe.hel.complete <- hclust(spe.hel, method = "complete")
        plot(spe.hel.complete, main = "Helinger - Complete linkage")
        spe.hel.comp.coph <- cophenetic(spe.hel.complete)
        cor(spe.hel, spe.hel.comp.coph)
        
        # Average method (BEST)
        spe.hel.UPGMA <- hclust(spe.hel, method = "average")
        plot(spe.hel.UPGMA, main = "Helinger - UPGMA")
        spe.hel.UPGMA.coph <- cophenetic(spe.hel.UPGMA)
        cor(spe.hel, spe.hel.UPGMA.coph)
        
        # Centroid method
        spe.hel.centroid <- hclust(spe.hel, method = "centroid")
        plot(spe.hel.centroid, main = "Helinger - Centroid")
        spe.hel.centroid.coph <- cophenetic(spe.hel.centroid)
        cor(spe.hel, spe.hel.centroid.coph)
        
        # Ward method
        spe.hel.ward <- hclust(spe.hel, method = "ward.D2")
        plot(spe.hel.ward,  main = "Helinger - Ward")
        spe.hel.ward.coph <- cophenetic(spe.hel.ward)
        cor(spe.hel, spe.hel.ward.coph)
        
        

library(NbClust)


Nb.complete <-NbClust(spe[,2:58], diss=spe.hel, distance = NULL, min.nc=1, max.nc=16, 
                   method = "average", index = "ch")

Nb.complete
plot(Nb.complete$All.index, xlab="number of clusters", ylab="Calinski and Harabasz index")

UPGMA.dend <- as.dendrogram(spe.hel.UPGMA)
plot(UPGMA.dend)

library(dendextend)

Nb.UPGMA<-NbClust(spe[,2:58], diss=spe.hel, distance = NULL, min.nc=2, max.nc=16, 
                  method = "average", index="ch")
Nb.UPGMA
plot(Nb.UPGMA$All.index, xlab ="number of clusters", ylab = "Calinski and Harabs index")

#convert to dendrogram
UPGMA.dend <- as.dendrogram(spe.hel.UPGMA)


#define colors and sort according to tips in dendrogram
colors_to_use <- Nb.UPGMA$Best.partition
colors_to_use<-colors_to_use[order.dendrogram(UPGMA.dend)]



#change color of tips
labels_colors(UPGMA.dend) <- colors_to_use
plot(UPGMA.dend)


#change color of branches
labels_colors(UPGMA.dend)<-1
UPGMA.dend <- UPGMA.dend %>% color_branches(k = 5)
plot(UPGMA.dend)

cla.hel.UPGMA.coph <- cophenetic(cla.hel.UPGMA)
cor(cla.hel, cla.hel.UPGMA.coph)

plot(cla.hel, cla.hel.ward.coph,
     xlab = "Chord distance",
     ylab = "Cophenetic distance",
     asp = 1, xlim = c(0, sqrt(2)),
     ylim = c(0, sqrt(2)),
     main = c("Single linkage", paste("Cophenetic correlation =", round(cor(cla.hel, cla.hel.ward.coph), 3))))
abline(0, 1)
lines(lowess(cla.hel, cla.hel.ward.coph), col = "red", lwd=3)

#Unconstrained ordination
    #CA
      
#Constrained ordination 
    # CCA
        
            
