library(tidyverse)


data.long <- bigdata %>%
  pivot_longer(!Animal_ID, names_to = "Variables", values_to = "Value") %>%
  as.data.frame()

data.long <- inner_join(data.long, bigmeta, by = "Animal_ID")

#Collect all outcome variables
data_names <- unique(data.long$Variables) 

#Create an empty container list of the appropriate length
data.list <- vector(mode = "list", length = length(data_names)) 

for(n in 1:length(data_names)){
  data.list[[n]] <- data.long %>% filter(Variables == data_names[n])
}
