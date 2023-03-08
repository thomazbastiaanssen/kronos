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
head(res_kronos)
```

    ##               p.val        r.sq         avg       acro  amplitude
    ## Species_1 0.7582688 0.019571468 -0.50018010  0.8700673 0.12469734
    ## Species_2 0.2466698 0.095143490  0.08157565 23.1023980 0.30710811
    ## Species_3 0.8405307 0.012332033  0.11704177 17.2228779 0.09742959
    ## Species_4 0.4882364 0.049921918  1.70138007 17.0889161 0.11382561
    ## Species_5 0.7613635 0.019286195  0.64496662 18.7856844 0.09737362
    ## Species_6 0.8692770 0.009956774 -0.09725603 22.8010299 0.09430708

``` r
head(res_jtk)
```

    ##           BH.Q     ADJ.P PER LAG        AMP
    ## Species_1    1 1.0000000  24  NA 0.23877198
    ## Species_2    1 0.4066034  24  21 0.23017044
    ## Species_3    1 1.0000000  24  NA 0.02678880
    ## Species_4    1 0.7528835  24  15 0.13805803
    ## Species_5    1 1.0000000  24  NA 0.01416228
    ## Species_6    1 1.0000000  24  NA 0.01824839

``` r
head(res_limo)
```

    ##      gene_id    time_cos    time_sin     AveExpr         F   P.Value adj.P.Val
    ## 1: Species_1  0.12147634  0.02815896 -0.49954083 0.3067721 0.7366005 0.9999392
    ## 2: Species_2  0.29866764 -0.07150548  0.07876214 1.8694925 0.1602636 0.9999392
    ## 3: Species_3 -0.01968560 -0.09542013  0.11410717 0.1931209 0.8247301 0.9999392
    ## 4: Species_4 -0.02689312 -0.11060303  1.69798652 0.4265757 0.6540805 0.9999392
    ## 5: Species_5  0.01988801 -0.09532097  0.64195753 0.2427789 0.7849677 0.9999392
    ## 6: Species_6  0.08969922 -0.02911831 -0.09833916 0.1732639 0.8412013 0.9999392

``` r
head(res_cosinor2)
```

    ##                    r Percent.rhythm    p.value  acrophase
    ## Species_1 0.13989806    0.019571468 0.45288540 -0.2277831
    ## Species_2 0.30845338    0.095143490 0.09135956 -6.0481937
    ## Species_3 0.11104969    0.012332033 0.55202111 -4.5089389
    ## Species_4 0.22343213    0.049921918 0.22695396 -4.4738678
    ## Species_5 0.13887474    0.019286195 0.45623024 -4.9180807
    ## Species_6 0.09978364    0.009956774 0.59328826 -5.9692957

``` r
res_tot = do.call(cbind, list(res_kronos, res_jtk, res_limo, res_cosinor2))
colnames(res_tot) = paste(colnames(res_tot), c(rep("kronos", ncol(res_kronos)),
                                               rep("jtk",    ncol(res_jtk)),
                                               rep("limo",    ncol(res_limo)),
                                               rep("cosinor2",    ncol(res_cosinor2))), sep = "_")
head(res_tot)
```

    ##           p.val_kronos r.sq_kronos  avg_kronos acro_kronos amplitude_kronos
    ## Species_1    0.7582688 0.019571468 -0.50018010   0.8700673       0.12469734
    ## Species_2    0.2466698 0.095143490  0.08157565  23.1023980       0.30710811
    ## Species_3    0.8405307 0.012332033  0.11704177  17.2228779       0.09742959
    ## Species_4    0.4882364 0.049921918  1.70138007  17.0889161       0.11382561
    ## Species_5    0.7613635 0.019286195  0.64496662  18.7856844       0.09737362
    ## Species_6    0.8692770 0.009956774 -0.09725603  22.8010299       0.09430708
    ##           BH.Q_jtk ADJ.P_jtk PER_jtk LAG_jtk    AMP_jtk gene_id_limo
    ## Species_1        1 1.0000000      24      NA 0.23877198    Species_1
    ## Species_2        1 0.4066034      24      21 0.23017044    Species_2
    ## Species_3        1 1.0000000      24      NA 0.02678880    Species_3
    ## Species_4        1 0.7528835      24      15 0.13805803    Species_4
    ## Species_5        1 1.0000000      24      NA 0.01416228    Species_5
    ## Species_6        1 1.0000000      24      NA 0.01824839    Species_6
    ##           time_cos_limo time_sin_limo AveExpr_limo    F_limo P.Value_limo
    ## Species_1    0.12147634    0.02815896  -0.49954083 0.3067721    0.7366005
    ## Species_2    0.29866764   -0.07150548   0.07876214 1.8694925    0.1602636
    ## Species_3   -0.01968560   -0.09542013   0.11410717 0.1931209    0.8247301
    ## Species_4   -0.02689312   -0.11060303   1.69798652 0.4265757    0.6540805
    ## Species_5    0.01988801   -0.09532097   0.64195753 0.2427789    0.7849677
    ## Species_6    0.08969922   -0.02911831  -0.09833916 0.1732639    0.8412013
    ##           adj.P.Val_limo r_cosinor2 Percent.rhythm_cosinor2 p.value_cosinor2
    ## Species_1      0.9999392 0.13989806             0.019571468       0.45288540
    ## Species_2      0.9999392 0.30845338             0.095143490       0.09135956
    ## Species_3      0.9999392 0.11104969             0.012332033       0.55202111
    ## Species_4      0.9999392 0.22343213             0.049921918       0.22695396
    ## Species_5      0.9999392 0.13887474             0.019286195       0.45623024
    ## Species_6      0.9999392 0.09978364             0.009956774       0.59328826
    ##           acrophase_cosinor2
    ## Species_1         -0.2277831
    ## Species_2         -6.0481937
    ## Species_3         -4.5089389
    ## Species_4         -4.4738678
    ## Species_5         -4.9180807
    ## Species_6         -5.9692957

``` r
plot(res_kronos$p.val, res_jtk$ADJ.P)
```

![](benchmarking_files/figure-gfm/assess%20tools-1.png)<!-- -->

``` r
hist(res_kronos$p.val, breaks = 20)
```

![](benchmarking_files/figure-gfm/assess%20tools-2.png)<!-- -->

``` r
hist(res_jtk$ADJ.P, breaks = 20)
```

![](benchmarking_files/figure-gfm/assess%20tools-3.png)<!-- -->

``` r
hist(res_cosinor2$p, breaks = 20)
```

![](benchmarking_files/figure-gfm/assess%20tools-4.png)<!-- -->

``` r
plot(res_kronos$amplitude, res_jtk$AMP)
```

![](benchmarking_files/figure-gfm/assess%20tools-5.png)<!-- -->

``` r
plot(res_kronos$acro, res_jtk$LAG)
```

![](benchmarking_files/figure-gfm/assess%20tools-6.png)<!-- -->

``` r
plot(res_kronos$acro, res_cosinor2$acrophase + 24)
```

![](benchmarking_files/figure-gfm/assess%20tools-7.png)<!-- -->
