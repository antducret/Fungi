# -*- coding: utf-8 -*-
# PREPARE SESSION ---------------------------------------------------------

# +
# Set directory ?
# Librairies ?
#setwd(dir = "/home/filsdufrere/Documents/Fungi/Projet_Rstudio")
#setwd("/home/anthoney/Documents/Master/MiR/Fungi/Projet_Rstudio/")
# -

library(tidyverse)
# library(vegan)

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
    colnames(env) <- c("Latitude", "Longitude", "Depth", "Sediment" , "O2" , "P04", "N03", "Temperature", "Salinity" , "SiO4")    # Rename ENV variables

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

# Remove duplicate sites and sums corresponding classes

    env = env[-c(10,35,39,52,53),]   
    cla[9,] = cla[9,]+cla[10,]
    cla[34,]= cla[34,]+cla[35,]
    cla[38,]= cla[38,] + cla[39,]
    cla[51,]= cla[51,] + cla[52,] + cla[53,]
    
    cla = cla[-c(10,35,39,52,53),]

# Reduce size of class name 

    colnames(cla) = c("Agarico.","Agaricostilbo.","Chytridio.","Dothideo.","Eurotio.","Glomero.","Lecanoro.","Leotio.","Microbotryo.","Monoblepharido.","Ochroconis","Pezizo.","Saccharo.","Sordario.","Taphrino.","Ustilagino.","Zygo.")

# Data exploration SPE  ---------------------------------------------------------------

#Summary
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
          geom_point(data = coo, aes(x=Longitude, y=Latitude, size = sites.pres,colour = "red"), inherit.aes = FALSE) +
          xlab("° Longitude") + ylab("° Latitude") +
          labs(colour ='Sites of study',size = "Number of different present classes" , title = "Location of sites")
        map

#Unsupervised classification -----------------------------------------------------------------------------------------------------------
  ## Hierarchical agglomerative clustering of the classes abundance
        library(vegan)
        # Compute Helinger distance
        cla <-  t(cla)
        cla.norm <- decostand(cla, "normalize")
        cla.hel <- vegdist(cla.norm, "hel")
  
        # Average method (BEST)
        
            cla.hel.UPGMA <- hclust(cla.hel, method = "average")
            plot(cla.hel.UPGMA, main = "Helinger - UPGMA")
            cla.hel.UPGMA.coph <- cophenetic(cla.hel.UPGMA)
            cor(cla.hel, cla.hel.UPGMA.coph)


  library(NbClust)
  par(mfrow = c(1, 1))
  library(dendextend)
  
  Nb.UPGMA<-NbClust(cla, diss=cla.hel, distance = NULL, min.nc=1, max.nc=16,
                    method = "average", index="ch")
  Nb.UPGMA
  
  #plot(Nb.UPGMA$All.index, xlab ="number of clusters", ylab = "Calinski and Harabs index")
  
  #convert to dendrogram
  plot(cla.hel.UPGMA, main = "Helinger - UPGMA")

  UPGMA.dend <- as.dendrogram(cla.hel.UPGMA)
  plot(UPGMA.dend)
  par(mfrow = c(1,1))
  #define colors and labels  and sort according to tips in dendrogram
  colors_to_use <- Nb.UPGMA$Best.partition
  colors_to_use<-colors_to_use[order.dendrogram(UPGMA.dend)]
  colors_to_use[colors_to_use!= 8] = seq(1,7)
  labels_to_use <- strtrim(cla.hel.UPGMA$labels, 11)
  labels_to_use <- labels_to_use[order.dendrogram(UPGMA.dend)]
  
  
  #change color of tips
  labels_colors(UPGMA.dend) <- colors_to_use
  labels(UPGMA.dend)<- labels_to_use
  UPGMA.dend <-color_branches(UPGMA.dend, col =  colors_to_use)
  plot(UPGMA.dend, main = "Helinger - UPGMA")


cla.hel.UPGMA.coph <- cophenetic(cla.hel.UPGMA)
cor(cla.hel, cla.hel.UPGMA.coph)

plot(cla.hel, cla.hel.UPGMA.coph,
     xlab = "Hellinger distance",
     ylab = "Cophenetic distance",
     asp = 1, xlim = c(0, sqrt(2)),
     ylim = c(0, sqrt(2)),
     main = c("Average linkage", paste("Cophenetic correlation =", round(cor(cla.hel, cla.hel.UPGMA.coph), 2))))
abline(0, 1)
lines(lowess(cla.hel, cla.hel.UPGMA.coph), col = "red", lwd=3)

?fviz_ca_biplot

columns

#CA (Correspondence analysis)
library(factoextra)
library("FactoMineR")
cla.ca2 <- CA(cla, graph = FALSE)
columns <- get_ca_col(cla.ca2)
p1 <- fviz_ca_biplot(cla.ca2, repel = TRUE, col.col = env$Sediment, label = 'row') +
scale_color_manual(name = "Sites", labels = c("Pelagic","Benthic"),values= c('orange','green'))
jpeg("CA_plot.jpg")
print(p1)
dev.off()
p1

?fviz_ca_biplot

#Constrained ordination ------------------------------------------------------------------------------------------------------------
## CCA of untransformed fish species data, constrained by all environmental variables (env3)
cla = t(cla)
cla.ca <- cca(cla, env)
summary(cla.ca)		# default scaling 2
summary(cla.ca, scaling = 1)

# Scree plot and broken stick model using vegan's screeplot.cca()
screeplot(cla.ca, npcs = length(cla.ca$CA$eig))

par(mfrow = c(1,1))
# Scaling 3: Compromise between both before TO COMMENT
plot(cla.ca, 
     scaling = 3,
     main = "CA fungi abundances - biplot scaling 3"
)


# ordisurf(cla.ca, env$sal, add = TRUE)
# ordisurf(cla.ca, env$sil, add = TRUE, col = "green")
