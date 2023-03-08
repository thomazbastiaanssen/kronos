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

    ##               p.val        r.sq         avg        acro amplitude
    ## Species_1 0.3820247 0.066424615 -0.49595787 22.22889539 0.2283776
    ## Species_2 0.2400488 0.096900343  0.07935109  0.00786307 0.3106057
    ## Species_3 0.8303727 0.013189438  0.11563778 22.22962151 0.0975930
    ## Species_4 0.4438078 0.056374569  1.69947939 22.63222907 0.1170285
    ## Species_5 0.0825625 0.163189071  0.63347595  5.05773263 0.2852202
    ## Species_6 0.9043273 0.007157404 -0.10048064  8.38252352 0.0810046

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
    ## 1: Species_1  0.12147634  0.02815896 -0.49954083 0.2952323 0.7449227 0.9955678
    ## 2: Species_2  0.29866764 -0.07150548  0.07876214 1.8959796 0.1548487 0.9955678
    ## 3: Species_3 -0.01968560 -0.09542013  0.11410717 0.1933117 0.8244931 0.9955678
    ## 4: Species_4 -0.02689312 -0.11060303  1.69798652 0.4814756 0.6191143 0.9955678
    ## 5: Species_5  0.01988801 -0.09532097  0.64195753 0.2603736 0.7712180 0.9955678
    ## 6: Species_6  0.08969922 -0.02911831 -0.09833916 0.1707561 0.8432411 0.9955678

``` r
head(res_cosinor2)
```

    ##                   F df1 df2         p    acrophase
    ## Species_1 0.9961109   2  28 0.3820247 -5.819511204
    ## Species_2 1.5021651   2  28 0.2400488 -0.002058547
    ## Species_3 0.1871201   2  28 0.8303727 -5.819701301
    ## Species_4 0.8363954   2  28 0.4438078 -5.925103715
    ## Species_5 2.7301830   2  28 0.0825625 -1.324111306
    ## Species_6 0.1009260   2  28 0.9043273 -2.194539525

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
