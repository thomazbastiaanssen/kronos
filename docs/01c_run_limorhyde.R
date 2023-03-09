#install.packages('limorhyde')
library(limorhyde)
library('annotate')
library('data.table')
library('limma')
library('limorhyde')


period = 24
qvalRhyCutoff = 0.15
qvalDrCutoff = 0.1

metadata_limo = cbind(metadata, limorhyde(metadata$Timepoint, 'time_'))

y = species
# rownames(y) <- y[,1]
# y <- y[,-1]
y <- y[,metadata_limo$woltka_ID]


res_limo = {
  design = model.matrix(~ time_cos + time_sin, data = metadata_limo)
  fit = lmFit(y, design)
  fit = eBayes(fit, trend = TRUE)
  rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
  setnames(rhyNow, 'rn', 'gene_id')
  rhyNow}

res_limo = res_limo[order(match(res_limo$gene_id,row.names(y))),]

colnames(res_limo)
colnames(res_limo) = c( "gene_id", "time_cos" , "time_sin",   "avg", "acro", "amplitude", "p.val")

rm(list=setdiff(ls(),  c("species", "metadata","res_jtk", "res_kronos", "res_limo", "res_cosinor2")))
