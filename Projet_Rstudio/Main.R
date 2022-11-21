
# PREPARE SESSION ---------------------------------------------------------
# Set directory ?
# Librairies ?
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
    spe <- spe[,-(0:11)]
    site_names <- colnames(spe[,-1], prefix = "")
    colnames(spe) <- append(1:56, "Tot", 0) # Pourquoi ?
  
    # Pelagic
    spe.pel <- read.csv(PATH_SPE_PELAGIC, sep = ',', quote = "", row.names = NULL); 
    spe.pel <- spe.pel[,-(0:11)];
    colnames(spe.pel) <- append(1:56, "Tot", 0);
  
    # Benthic
    spe.ben <- read.csv(PATH_SPE_BENTHIC, sep = ',', quote = "", row.names = NULL) 
    spe.ben <- spe.ben[,-(0:11)]
    colnames(spe.ben) <- append(1:56, "Tot", 0)
  
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
         labs(colour ='Sites of study',size = "Number of occurrences" , title = "Location of sites")
map 
    
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
          spe.pres <- apply(spe > 0, 1, sum)
          sites.pres <- apply(spe > 0, 2, sum)
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
         

# Plot species on map ?? 
            
