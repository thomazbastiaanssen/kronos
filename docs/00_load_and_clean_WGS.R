library(deleuze)

# species = read.delim("docs/wgs_species.csv")
# metadata = read.delim("docs/wgs_metadata.txt", sep = ",", row.names = 1)

species = read.delim("wgs_species.csv", row.names = 1)
metadata = read.delim("wgs_metadata.txt", sep = ",", row.names = 1)

species   <- species[apply(species == 0, 1, sum) <= round(ncol(species) * 0.9), ]   #remove rows with 10% or fewer hits
species   <- data.frame(genefilter::varFilter(as.matrix(species), var.cutoff = 0.4))

#CLR transform the species count table.
species <- getTableMeans(species) 
species <- species[,metadata$woltka_ID]


