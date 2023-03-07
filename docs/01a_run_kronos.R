library(tidyverse)

row.names(species) = species[, 1]
species            = species[,-1]

out_list = fw_kronos(species, formula = ~ time(Timepoint), metadata = metadata, pairwise = F)

res_kronos = do.call(rbind, lapply(FUN = getKronos_groupwise, X = out_list))[,-1]


rm(list=setdiff(ls(), c("species", "metadata","res_jtk", "res_kronos")))

