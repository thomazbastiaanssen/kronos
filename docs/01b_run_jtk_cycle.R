source("https://raw.githubusercontent.com/mfcovington/jtk-cycle/develop/JTK_CYCLE.R")

project <- "benchmark"

data = species
# rownames(data) <- data[,1]
# data <- data[,-1]
data <- data[,metadata$woltka_ID]
jtkdist(4, c(8, 8, 7, 8))

periods <- 4
jtk.init(periods,6)



  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  res_jtk = res
  
  colnames(res_jtk) 
  colnames(res_jtk) = c(  "q.val", "p.val" , "period",   "acro", "amplitude")
  rm(list=setdiff(ls(),  c("species", "metadata","res_jtk", "res_kronos", "res_limo", "res_cosinor2")))
  
