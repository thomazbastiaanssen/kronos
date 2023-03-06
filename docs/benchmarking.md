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
#devtools::install_github("thomazbastiaanssen/kronos") 
library(kronos)

#install.packages('limorhyde')
library(limorhyde)
```

``` r
source("00_load_and_clean_WGS.R")
```

``` r
source("01a_run_kronos.R")
source("01b_run_jtk_cycle.R")
source("01c_run_limorhyde.R")
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: XML

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     shift

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, second

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

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

    ##           BH.Q     ADJ.P PER LAG       AMP
    ## Species_1    1 0.7062927  24  21 0.3690645
    ## Species_2    1 0.5045277  24  21 0.3009061
    ## Species_3    1 1.0000000  24  NA 0.1064677
    ## Species_4    1 0.2168050  24  18 0.1665836
    ## Species_5    1 0.1992868  24   0 0.1001631
    ## Species_6    1 1.0000000  24  NA 0.0365927

``` r
head(res_limo)
```

    ##        gene_id   time_cos   time_sin    AveExpr         F      P.Value
    ## 1: Species_242 -0.4854064  0.9645182  1.4533285 33.056781 1.962474e-12
    ## 2: Species_286 -0.3284821  0.7049102  3.9795961 20.995880 1.121822e-08
    ## 3: Species_415 -0.1662919  0.4487244  1.8714063  9.141760 1.879773e-04
    ## 4: Species_263 -0.2576868 -0.4778930 -1.4230559  8.975189 2.177556e-04
    ## 5: Species_452 -0.2576868 -0.4778930 -1.4230559  8.975189 2.177556e-04
    ## 6: Species_336 -0.5691492 -0.1913979 -0.4585398  7.586936 7.510226e-04
    ##       adj.P.Val
    ## 1: 9.812371e-10
    ## 2: 2.804555e-06
    ## 3: 2.177556e-02
    ## 4: 2.177556e-02
    ## 5: 2.177556e-02
    ## 6: 5.644168e-02

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
