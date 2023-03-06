``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0      ✔ purrr   1.0.1 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.5.0 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(kronos)
data("kronos_demo")
```

``` r
source("00_load_and_clean_WGS.R")
```

``` r
source("01a_run_kronos.R")
source("01b_run_jtk_cycle.R")

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

    ##           BH.Q     ADJ.P PER LAG       AMP
    ## Species_1    1 0.7062927  24  21 0.3690645
    ## Species_2    1 0.5045277  24  21 0.3009061
    ## Species_3    1 1.0000000  24  NA 0.1064677
    ## Species_4    1 0.2168050  24  18 0.1665836
    ## Species_5    1 0.1992868  24   0 0.1001631
    ## Species_6    1 1.0000000  24  NA 0.0365927

``` r
plot(res_kronos$p.val, res_jtk$ADJ.P)
```

![](benchmarking_files/figure-gfm/run%20tools-1.png)<!-- -->

``` r
hist(res_kronos$p.val, breaks = 20)
```

![](benchmarking_files/figure-gfm/run%20tools-2.png)<!-- -->

``` r
hist(res_jtk$ADJ.P, breaks = 20)
```

![](benchmarking_files/figure-gfm/run%20tools-3.png)<!-- -->

``` r
plot(res_kronos$amplitude, res_jtk$AMP)
```

![](benchmarking_files/figure-gfm/run%20tools-4.png)<!-- -->

``` r
plot(res_kronos$acro, res_jtk$LAG)
```

![](benchmarking_files/figure-gfm/run%20tools-5.png)<!-- -->
