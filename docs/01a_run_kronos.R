library(tidyverse)

data.long <- species %>%
  pivot_longer(!Name, names_to = "woltka_ID", values_to = "Value") %>%
  as.data.frame()

data.long <- inner_join(data.long, metadata, by = "woltka_ID")

#Collect all outcome variables
data_names <- unique(data.long$Name) 

#Create an empty container list of the appropriate length
data.list <- vector(mode = "list", length = length(data_names)) 

for(n in 1:length(data_names)){
  data.list[[n]] <- data.long %>% filter(Name == data_names[n])
}

out_list = lapply(X = data.list, 
                  FUN = function(y){
                    kronos(data = y, 
                           Value ~ time(Timepoint), 
                           period = 24, 
                           pairwise = F, verbose = F)
                  }
)
names(out_list) <- data_names

res_kronos = do.call(rbind, lapply(FUN = getKronos_groupwise, X = out_list))[,-1]

rm(list=setdiff(ls(), c("species", "metadata","res_jtk", "res_kronos")))

