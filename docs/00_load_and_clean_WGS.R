library(deleuze)

species = read.delim("wgs_species.csv")
metadata = read.delim("wgs_metadata.txt", sep = ",", row.names = 1)

#CLR transform the species count table.
species[,-1] <-  getTableMeans(species[,-1]) 
species[,-1] <- species[,-1][,metadata$woltka_ID]

species = species[1:500,]

