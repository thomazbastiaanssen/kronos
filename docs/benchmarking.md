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
```

    ## [1] "BH.Q"  "ADJ.P" "PER"   "LAG"   "AMP"

``` r
source("01c_run_limorhyde.R")
```

    ## [1] "gene_id"   "time_cos"  "time_sin"  "AveExpr"   "F"         "P.Value"  
    ## [7] "adj.P.Val"

``` r
source("01d_run_cosinor2.R")
```

    ## [1] "r"                                                                          
    ## [2] "Percent.rhythm"                                                             
    ## [3] "p.value"                                                                    
    ## [4] "do.call.rbind..lapply.X...fw_cosinor.x...data..formula....time.Timepoint..."

``` r
head(res_kronos)
```

    ##               p.val       r.sq         avg      acro amplitude
    ## Species_2 0.3515879 0.07194470 -0.59007048 22.658650 0.2627859
    ## Species_3 0.7134866 0.02382527 -0.55378379 18.658285 0.1366910
    ## Species_4 0.2495125 0.09440261  1.03005782 19.273910 0.1397519
    ## Species_5 0.3632017 0.06978789 -0.03594562  5.880808 0.1856407
    ## Species_6 0.7985742 0.01593787 -0.76990221 12.854223 0.1140083
    ## Species_8 0.4024962 0.06293720  0.05263352 10.507180 0.1747327

``` r
head(res_jtk)
```

    ##               q.val     p.val period acro  amplitude
    ## Species_2 1.0000000 0.3774844     24   18 0.26479047
    ## Species_3 1.0000000 1.0000000     24   NA 0.04940765
    ## Species_4 0.9651211 0.1992868     24   15 0.13129342
    ## Species_5 1.0000000 1.0000000     24   NA 0.12170528
    ## Species_6 1.0000000 1.0000000     24   NA 0.10195058
    ## Species_8 1.0000000 1.0000000     24   NA 0.03600578

``` r
head(res_limo)
```

    ##      gene_id     time_cos    time_sin         avg         F     p.val     q.val
    ## 1: Species_2  0.246748747 -0.09039613 -0.59337082 1.2529959 0.2943813 0.8330279
    ## 2: Species_3  0.023440712 -0.13466616 -0.55802579 0.3689815 0.6932880 0.9590520
    ## 3: Species_4  0.045749257 -0.13205147  1.02585356 0.9832489 0.3811448 0.8767841
    ## 4: Species_5  0.005791869  0.18555031 -0.03017543 1.0099758 0.3714700 0.8707995
    ## 5: Species_6 -0.111169227 -0.02528427 -0.77047212 0.2338499 0.7923334 0.9770689
    ## 6: Species_8 -0.161557335  0.06656371  0.05502427 0.9355144 0.3990815 0.8796326

``` r
head(res_cosinor2)
```

    ##                   r       r.sq      p.val      acro
    ## Species_2 0.2682251 0.07194470 0.14458386 -5.932021
    ## Species_3 0.1543544 0.02382527 0.40706041 -4.884728
    ## Species_4 0.3072501 0.09440261 0.09269616 -5.045898
    ## Species_5 0.2641740 0.06978789 0.15098209 -1.539592
    ## Species_6 0.1262453 0.01593787 0.49857470 -3.365228
    ## Species_8 0.2508729 0.06293720 0.17342296 -2.750773

``` r
res_tot = do.call(cbind, list(res_kronos, res_jtk, res_limo, res_cosinor2))
colnames(res_tot) = paste(colnames(res_tot), c(rep("kronos", ncol(res_kronos)),
                                               rep("jtk",    ncol(res_jtk)),
                                               rep("limo",    ncol(res_limo)),
                                               rep("cosinor2",    ncol(res_cosinor2))), sep = "_")

colnames(res_tot)
```

    ##  [1] "p.val_kronos"     "r.sq_kronos"      "avg_kronos"       "acro_kronos"     
    ##  [5] "amplitude_kronos" "q.val_jtk"        "p.val_jtk"        "period_jtk"      
    ##  [9] "acro_jtk"         "amplitude_jtk"    "gene_id_limo"     "time_cos_limo"   
    ## [13] "time_sin_limo"    "avg_limo"         "F_limo"           "p.val_limo"      
    ## [17] "q.val_limo"       "r_cosinor2"       "r.sq_cosinor2"    "p.val_cosinor2"  
    ## [21] "acro_cosinor2"

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
  # pivot_longer(!ID) %>% 
  # separate(name, into = c("parameter", "program"), sep = "_")  %>% 
  # filter(parameter %in% c("acro", "p.val")) %>% 
  # pivot_wider(names_from = parameter, values_from = value) %>% 
  
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
  UpSetR::upset(., order.by = "freq")
```

![](benchmarking_files/figure-gfm/plot%20tools-3.png)<!-- -->
