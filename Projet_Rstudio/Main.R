
# PATH ---------------------------------------------------------------
PATH_OCCURRENCE_TOT = "../Data/occurence_tot.csv"
PATH_OCCURRENCE_BENTHIC = "../Data/occurence_benthic.csv"
PATH_OCCURRENCE_PELAGIC = "../Data/occurence_pelagic.csv"
PATH_ENV = "../Data/env.csv"
# Load data ---------------------------------------------------------------

#Environmental parameters
env <- read.csv(PATH_ENV, sep = ',', skip = 1 )
env <- env[2205:2214]
env
#Occurrences
occ.tot <- read.csv(PATH_OCCURRENCE_TOT, sep = ',', row.names = NULL)
name_species <- occ.tot[,2:11]
occ.tot <- occ.tot[,-(0:11)]
site_names <- colnames(occ.tot[,-1], prefix = "")
colnames(occ.tot) <- append(1:56, "Tot", 0)


occ.pel <- read.csv(PATH_OCCURRENCE_PELAGIC, sep = ',', quote = "", row.names = NULL) 
occ.pel <- occ.pel[,-(0:11)]
colnames(occ.pel) <- append(1:56, "Tot", 0)

occ.ben <- read.csv(PATH_OCCURRENCE_BENTHIC, sep = ',', quote = "", row.names = NULL) 
occ.ben <- occ.ben[,-(0:11)]
colnames(occ.ben) <- append(1:56, "Tot", 0)


env
occ.tot
occ.pel 
occ.ben

