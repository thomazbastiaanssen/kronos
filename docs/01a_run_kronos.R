library(tidyverse)
library(kronos)

data = species
row.names(data) = data[, 1]
data            = data[,-1]

out_list = fw_kronos(data, formula = ~ time(Timepoint), metadata = metadata, pairwise = F)

res_kronos = do.call(rbind, lapply(FUN = getKronos_groupwise, X = out_list))[,-1]


rm(list=setdiff(ls(),  c("species", "metadata","res_jtk", "res_kronos", "res_limo", "res_cosinor2")))

