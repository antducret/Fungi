
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
    env <- env[2205:2214]                                                                         # Keep only environmental parameters
    env$sediment <- as.factor(env$sediment)
    colnames(env) <- c("lat", "long", "depth", "sed" , "disO2" , "P", "N", "T", "sal" , "sil")    # Rename ENV variables

  # SPE
    spe <- read.csv(PATH_SPE, sep = ',', row.names = NULL)
    name_species <- spe[,2:11]                                                                    # Save nomenclatures of all species

    spe <- spe[spe$Class != '',]                                                                  # Keep only species with specified class
    spe$Class <- as.factor(spe$Class)

    cla <- aggregate(spe[,-(0:12)], by = list(spe$Class), FUN = "sum")                            # Define CLA as SPE without nomenclature

    rownames(cla) <- cla[,1]
    cla <-  cla[,-1]
    colnames(cla) <- 1:56
    cla <-  t(cla)

  # remove "empty" sites from CLA and ENV
    # total presence
      sites.pres <- apply(cla> 0, 1, sum)

    # remove empty sites
      index = sites.pres > 0
      cla <- cla[index,]
      env <- env[index,]

# Data exploration SPE  ---------------------------------------------------------------

    # Summary
      summary(cla)
      str(cla)

    # relative presence of each classes

        # total presence for each classes
          cla.pres <- apply(cla> 0, 2, sum)
          sites.pres <- apply(cla> 0, 1, sum)

    # Proportion of zeros in the community data set
        sum(cla == 0) / (nrow(cla) * ncol(cla))

    # Plot barplots
        # Sort the results in increasing order
        sort(cla.pres)


        barplot(sites.pres,
               las = 1,
               xlab = "Sites",
               ylab = "sites richness",
               col = gray(5 : 0 / 5),
               horiz=F,
       )

        barplot(cla.pres,
                las = 1,
                xlab = "Sites",
                ylab = "classies richness",
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
          geom_point(data = coo, aes(x=long, y=lat, size = sites.pres,colour = "red"), inherit.aes = FALSE) +
          xlab("° Longitude") + ylab("° Latitude") +
          labs(colour ='Sites of study',size = "Number of occurrences of classes" , title = "Location of sites")
        map



#Unsupervised classification -----------------------------------------------------------------------------------------------------------
  ## Hierarchical agglomerative clustering of the classes abundance

        # Compute Helinger distance
        cla <-  t(cla)
        cla.norm <- decostand(cla, "normalize")
        cla.hel <- vegdist(cla.norm, "hel")

        # METHOD
          # CREATE DENDROGRAM
          # PLOT DENDROGRAM
          # COMPUTE COPHENETIC MATRIX
          # COMPUTE CORRELATION

        # Single method
        cla.hel.single <- hclust(cla.hel, method = "single")
        plot(cla.hel.single, main = "Helinger - Single linkage")
        cla.hel.single.coph <- cophenetic(cla.hel.single)
        cor(cla.hel, cla.hel.single.coph)

        # Complete method
        cla.hel.complete <- hclust(cla.hel, method = "complete")
        plot(cla.hel.complete, main = "Helinger - Complete linkage")
        cla.hel.comp.coph <- cophenetic(cla.hel.complete)
        cor(cla.hel, cla.hel.comp.coph)

        # Average method (BEST)
        cla.hel.UPGMA <- hclust(cla.hel, method = "average")
        plot(cla.hel.UPGMA, main = "Helinger - UPGMA")
        cla.hel.UPGMA.coph <- cophenetic(cla.hel.UPGMA)
        cor(cla.hel, cla.hel.UPGMA.coph)

        # Centroid method
        cla.hel.centroid <- hclust(cla.hel, method = "centroid")
        plot(cla.hel.centroid, main = "Helinger - Centroid")
        cla.hel.centroid.coph <- cophenetic(cla.hel.centroid)
        cor(cla.hel, cla.hel.centroid.coph)

        # Ward method
        cla.hel.ward <- hclust(cla.hel, method = "ward.D2")
        plot(cla.hel.ward,  main = "Helinger - Ward")
        cla.hel.ward.coph <- cophenetic(cla.hel.ward)
        cor(cla.hel, cla.hel.ward.coph)



library(NbClust)


Nb.complete <-NbClust(cla, diss=cla.hel, distance = NULL, min.nc=1, max.nc=16,
                   method = "average", index = "ch")

Nb.complete
plot(Nb.complete$All.index, xlab="number of clusters", ylab="Calinski and Harabasz index")

UPGMA.dend <- as.dendrogram(cla.hel.UPGMA)
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

cla <-t(cla)
cla.ca <- cca(cla)
summary(cla.ca)		# default scaling 2
summary(cla.ca, scaling = 1)

# Scree plot and broken stick model using vegan's screeplot.cca()
screeplot(cla.ca, bstick = TRUE, npcs = length(cla.ca$CA$eig))

# CA biplots
par(mfrow = c(1, 2))
# Scaling 1: sites are centroids of classes
plot(cla.ca,
     scaling = 1,
     main = "CA fungi abundances - biplot scaling 1"
)
# Scaling 2 (default): classes are centroids of sites
plot(cla.ca, main = "CA fungi abundances - biplot scaling 2")


# Curve fitting in a CA biplot
plot(cla.ca, main = "CA fungi abundances - scaling 2",
     sub = "Fitted curves: discharge (red), ammonium (green)")
cla.ca.env <- envfit(cla.ca ~ long+ lat + depth + T + P + N + sal + sil + disO2 , env)
plot(cla.ca.env)  # Two arrows·

#ordisurf(cla.ca, env$sal, add = TRUE)
#ordisurf(cla.ca, env$sil, add = TRUE, col = "green")


#Constrained ordination ------------------------------------------------------------------------------------------------------------
    # CCA
