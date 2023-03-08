library(tidyverse)
library(cosinor)
library(cosinor2)

data = species
row.names(data) = data[, 1]
data            = data[,-1]
data            = data[,metadata$woltka_ID]

fw_cosinor = function (x, formula, metadata, period = 24, verbose = F, 
          pairwise = F) 
{
  formula = update.formula(x_feature ~ ., formula)
  apply(X = x, MARGIN = 1, FUN = function(y) {
    cosinor:::cosinor.lm(formula = formula, data = cbind(x_feature = y, 
                                           metadata), period = period)
 }, simplify = F)
}

res_cosinor2 = data.frame(cbind(do.call(rbind, 
                   lapply(X = fw_cosinor(x = data, formula =  ~ time(Timepoint), period = 24, metadata = metadata), FUN = cosinor.PR)),
                 do.call(rbind, 
                         lapply(X = fw_cosinor(x = data, formula =  ~ time(Timepoint), period = 24, metadata = metadata), FUN = correct.acrophase))))

colnames(res_cosinor2)[4] = "acrophase"

rm(list=setdiff(ls(),  c("species", "metadata","res_jtk", "res_kronos", "res_limo", "res_cosinor2")))


