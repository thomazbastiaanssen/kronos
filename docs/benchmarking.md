``` r
library(tidyverse)

#devtools::install_github("thomazbastiaanssen/kronos") 
library(kronos)

#install.packages('limorhyde')
library(limorhyde)

library(cosinor)
library(cosinor2)
```

``` r
source("00_load_and_clean_WGS.R")
```

``` r
source("01a_run_kronos.R")
source("01b_run_jtk_cycle.R")
source("01c_run_limorhyde.R")
source("01d_run_cosinor2.R")
```

``` r
res_tot = do.call(cbind, list(res_kronos, res_jtk, res_limo, res_cosinor2))
colnames(res_tot) = paste(colnames(res_tot), c(rep("kronos", ncol(res_kronos)),
                                               rep("jtk",    ncol(res_jtk)),
                                               rep("limo",    ncol(res_limo)),
                                               rep("cosinor2",    ncol(res_cosinor2))), sep = "_")
```

``` r
par(mfcol = c(2, 2))

hist(res_kronos$p.val,   breaks = 20)
hist(res_jtk$p.val,      breaks = 20)
hist(res_cosinor2$p.val, breaks = 20)
hist(res_limo$p.val,     breaks = 20)
```

![](benchmarking_files/figure-gfm/plot%20stats-1.png)<!-- -->

``` r
res_tot %>%
  rownames_to_column("ID") %>%
  dplyr::select(!gene_id_limo) %>% 
  group_by(ID) %>% 
  filter(min(p.val_kronos,p.val_limo,p.val_jtk,p.val_cosinor2) < 0.05) %>% 
  ungroup() %>% 
  
  ggplot() +
  aes(x = .panel_x, y = .panel_y) +
  geom_point(aes(x = .panel_x, y = .panel_y)) +
  ggforce::facet_matrix(vars(acro_kronos, acro_jtk, acro_cosinor2)) +
  theme_bw() +
  ggtitle("Concordance in acrophase detection between Kronos, JTK_CYCLE and Cosinor2")
```

![](benchmarking_files/figure-gfm/plot%20tools-1.png)<!-- -->

``` r
list(kronos    = row.names(res_tot[res_tot$p.val_kronos   < 0.05,]),
     jtk_cycle = row.names(res_tot[res_tot$p.val_jtk      < 0.05,]),
     limorhyde = row.names(res_tot[res_tot$p.val_limo     < 0.05,]),
     cosinor2  = row.names(res_tot[res_tot$p.val_cosinor2 < 0.05,])) %>% 
  ggvenn::ggvenn()
```

![](benchmarking_files/figure-gfm/plot%20tools-2.png)<!-- -->

``` r
list(kronos    = row.names(res_tot[res_tot$p.val_kronos   < 0.05,]),
     jtk_cycle = row.names(res_tot[res_tot$p.val_jtk      < 0.05,]),
     limorhyde = row.names(res_tot[res_tot$p.val_limo     < 0.05,]),
     cosinor2  = row.names(res_tot[res_tot$p.val_cosinor2 < 0.05,])) %>% 
  UpSetR::fromList() %>% 
  UpSetR::upset(., order.by = "freq", 
                mainbar.y.label = "Number of p-values < 0.05\n per intersection", 
                sets.x.label = "Number of p-values < 0.05\n per program")
```

![](benchmarking_files/figure-gfm/plot%20tools-3.png)<!-- -->
