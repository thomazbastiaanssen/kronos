## Preparation

Here, we will compare the performance of Kronos to that of three common
R-based tools to assess rhythmic data. First, we will load the required
packages.

``` r
#Wrangling
library(tidyverse)

#Plotting
library(ggforce)
library(ggvenn)
library(UpSetR)

#Load circadian rhythm packages
library(kronos) #devtools::install_github("thomazbastiaanssen/kronos") 
library(limorhyde)
library(cosinor)
library(cosinor2)
```

Now we will load and prepare our high dimensional microbiome data for
assessment.

``` r
source("00_load_and_clean_WGS.R")
```

Next, we apply the four packages to the data. Code to run JTK_CYCLE,
Limorhyde and Cosinor2 was adapted from their respective vignettes.

``` r
source("01a_run_kronos.R")
source("01b_run_jtk_cycle.R")
source("01c_run_limorhyde.R")
source("01d_run_cosinor2.R")
```

Here, we collate all outcome parameters into a single object to compare
them.

``` r
res_tot = do.call(cbind, list(res_kronos, res_jtk, res_limo, res_cosinor2))
colnames(res_tot) = paste(colnames(res_tot), c(rep("kronos", ncol(res_kronos)),
                                               rep("jtk",    ncol(res_jtk)),
                                               rep("limo",    ncol(res_limo)),
                                               rep("cosinor2",    ncol(res_cosinor2))), sep = "_")
```

## Comparing p-value distributions between programs

Theoretically, we expect a uniform distribution between 0-1 (for the
cases where H_0 is true) as well as an inflation of low p-values (for
the cases where H_1 is true).

``` r
par(mfcol = c(2, 2))

hist(res_kronos$p.val,   breaks = 20)
hist(res_jtk$p.val,      breaks = 20)
hist(res_cosinor2$p.val, breaks = 20)
hist(res_limo$p.val,     breaks = 20)
mtext("P-value histograms", side = 3, line = -1, outer = TRUE)
```

![](benchmarking_files/figure-gfm/plot%20stats-1.png)<!-- -->

Both Limorhyde and Kronos show a desirable distribution of p-values,
whereas JTK_CYCLE has an inflation on the higher end and Cosinor2 has a
skewed, non-uniform distribution in the H_0 region.

## Comparing acrophases between programs

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

![](benchmarking_files/figure-gfm/plot%20acrophases-1.png)<!-- --> All
three programs that produce acrophases are in agreement as to where the
acrophases of features foudn to be rhythmic at p \< 0.05 at least once
are. Limorhyde does not have a clear method to return acrophases.

## Plot agreement in hits between programs

In the venn diagram we see how many features are considered hits (p\<
0.05) by the respective programs Alternatively, the UpSetR plot below
shows the same.

``` r
list(kronos    = row.names(res_tot[res_tot$p.val_kronos   < 0.05,]),
     jtk_cycle = row.names(res_tot[res_tot$p.val_jtk      < 0.05,]),
     limorhyde = row.names(res_tot[res_tot$p.val_limo     < 0.05,]),
     cosinor2  = row.names(res_tot[res_tot$p.val_cosinor2 < 0.05,])) %>% 
  ggvenn::ggvenn() + 
  ggtitle("Venn diagram showing agreement between features rhythmic at p < 0.05")
```

![](benchmarking_files/figure-gfm/plot%20hist-1.png)<!-- -->

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

![](benchmarking_files/figure-gfm/plot%20hist-2.png)<!-- --> From this,
we can conclude that Kronos performs as well as (or better than) the
best comparable program in each category. Furthermore, Kronos allows for
convenient testing of arbitrarily complex models and comes with
pre-designed ggplot-compatible plotting functions.

## Session Info

``` r
sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.1.2 (2021-11-01)
    ##  os       Ubuntu 20.04.3 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language en_GB:en
    ##  collate  en_GB.UTF-8
    ##  ctype    en_GB.UTF-8
    ##  tz       Europe/Dublin
    ##  date     2023-03-17
    ##  pandoc   2.19.2 @ /usr/lib/rstudio/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package          * version  date (UTC) lib source
    ##  annotate           1.72.0   2021-10-26 [1] Bioconductor
    ##  AnnotationDbi      1.56.2   2021-11-09 [1] Bioconductor
    ##  assertthat         0.2.1    2019-03-21 [1] CRAN (R 4.1.2)
    ##  backports          1.4.1    2021-12-13 [1] CRAN (R 4.1.2)
    ##  base64enc          0.1-3    2015-07-28 [1] CRAN (R 4.1.2)
    ##  Biobase            2.54.0   2021-10-26 [1] Bioconductor
    ##  BiocGenerics       0.40.0   2021-10-26 [1] Bioconductor
    ##  Biostrings         2.62.0   2021-10-26 [1] Bioconductor
    ##  bit                4.0.4    2020-08-04 [1] CRAN (R 4.1.2)
    ##  bit64              4.0.5    2020-08-30 [1] CRAN (R 4.1.2)
    ##  bitops             1.0-7    2021-04-24 [1] CRAN (R 4.1.2)
    ##  blob               1.2.3    2022-04-10 [1] CRAN (R 4.1.2)
    ##  broom              0.8.0    2022-04-13 [1] CRAN (R 4.1.2)
    ##  cachem             1.0.6    2021-08-19 [1] CRAN (R 4.1.2)
    ##  cellranger         1.1.0    2016-07-27 [1] CRAN (R 4.1.2)
    ##  checkmate          2.0.0    2020-02-06 [1] CRAN (R 4.1.2)
    ##  cli                3.3.0    2022-04-25 [1] CRAN (R 4.1.2)
    ##  cluster            2.1.2    2021-04-17 [4] CRAN (R 4.0.5)
    ##  colorspace         2.0-3    2022-02-21 [1] CRAN (R 4.1.2)
    ##  cosinor          * 1.2.2    2022-05-24 [1] CRAN (R 4.1.2)
    ##  cosinor2         * 0.2.1    2018-10-15 [1] CRAN (R 4.1.2)
    ##  cowplot            1.1.1    2020-12-30 [1] CRAN (R 4.1.2)
    ##  crayon             1.5.1    2022-03-26 [1] CRAN (R 4.1.2)
    ##  data.table         1.14.2   2021-09-27 [1] CRAN (R 4.1.2)
    ##  DBI                1.1.2    2021-12-20 [1] CRAN (R 4.1.2)
    ##  dbplyr             2.2.0    2022-06-05 [1] CRAN (R 4.1.2)
    ##  deleuze          * 0.1.0.0  2022-12-30 [1] Github (thomazbastiaanssen/deleuze@2c19a87)
    ##  digest             0.6.29   2021-12-01 [1] CRAN (R 4.1.2)
    ##  dplyr            * 1.0.9    2022-04-28 [1] CRAN (R 4.1.2)
    ##  ellipsis           0.3.2    2021-04-29 [1] CRAN (R 4.1.2)
    ##  evaluate           0.15     2022-02-18 [1] CRAN (R 4.1.2)
    ##  fansi              1.0.3    2022-03-24 [1] CRAN (R 4.1.2)
    ##  farver             2.1.1    2022-07-06 [1] CRAN (R 4.1.2)
    ##  fastmap            1.1.0    2021-01-25 [1] CRAN (R 4.1.2)
    ##  forcats          * 0.5.1    2021-01-27 [1] CRAN (R 4.1.2)
    ##  foreign            0.8-81   2020-12-22 [4] CRAN (R 4.0.3)
    ##  Formula            1.2-4    2020-10-16 [1] CRAN (R 4.1.2)
    ##  fs                 1.5.2    2021-12-08 [1] CRAN (R 4.1.2)
    ##  genefilter         1.76.0   2021-10-26 [1] Bioconductor
    ##  generics           0.1.2    2022-01-31 [1] CRAN (R 4.1.2)
    ##  GenomeInfoDb       1.30.1   2022-01-30 [1] Bioconductor
    ##  GenomeInfoDbData   1.2.7    2021-12-25 [1] Bioconductor
    ##  ggforce          * 0.3.3    2021-03-05 [1] CRAN (R 4.1.2)
    ##  ggplot2          * 3.3.6    2022-05-03 [1] CRAN (R 4.1.2)
    ##  ggvenn           * 0.1.9    2021-06-29 [1] CRAN (R 4.1.2)
    ##  glue               1.6.2    2022-02-24 [1] CRAN (R 4.1.2)
    ##  gridExtra          2.3      2017-09-09 [1] CRAN (R 4.1.2)
    ##  gtable             0.3.0    2019-03-25 [1] CRAN (R 4.1.2)
    ##  haven              2.5.0    2022-04-15 [1] CRAN (R 4.1.2)
    ##  highr              0.9      2021-04-16 [1] CRAN (R 4.1.2)
    ##  Hmisc              4.6-0    2021-10-07 [1] CRAN (R 4.1.2)
    ##  hms                1.1.1    2021-09-26 [1] CRAN (R 4.1.2)
    ##  htmlTable          2.4.0    2022-01-04 [1] CRAN (R 4.1.2)
    ##  htmltools          0.5.2    2021-08-25 [1] CRAN (R 4.1.2)
    ##  htmlwidgets        1.5.4    2021-09-08 [1] CRAN (R 4.1.2)
    ##  httpuv             1.6.5    2022-01-05 [1] CRAN (R 4.1.2)
    ##  httr               1.4.3    2022-05-04 [1] CRAN (R 4.1.2)
    ##  IRanges            2.28.0   2021-10-26 [1] Bioconductor
    ##  jpeg               0.1-9    2021-07-24 [1] CRAN (R 4.1.2)
    ##  jsonlite           1.8.0    2022-02-22 [1] CRAN (R 4.1.2)
    ##  KEGGREST           1.34.0   2021-10-26 [1] Bioconductor
    ##  knitr              1.39     2022-04-26 [1] CRAN (R 4.1.2)
    ##  kronos           * 0.1.0.0  2023-03-07 [1] Github (thomazbastiaanssen/kronos@540e437)
    ##  labeling           0.4.2    2020-10-20 [1] CRAN (R 4.1.2)
    ##  later              1.3.0    2021-08-18 [1] CRAN (R 4.1.2)
    ##  lattice            0.20-45  2021-09-22 [4] CRAN (R 4.1.1)
    ##  latticeExtra       0.6-29   2019-12-19 [1] CRAN (R 4.1.2)
    ##  lifecycle          1.0.1    2021-09-24 [1] CRAN (R 4.1.2)
    ##  limorhyde        * 1.0.1    2022-02-18 [1] CRAN (R 4.1.2)
    ##  lubridate          1.8.0    2021-10-07 [1] CRAN (R 4.1.2)
    ##  magrittr           2.0.3    2022-03-30 [1] CRAN (R 4.1.2)
    ##  MASS               7.3-54   2021-05-03 [4] CRAN (R 4.0.5)
    ##  Matrix             1.4-0    2021-12-08 [4] CRAN (R 4.1.2)
    ##  matrixStats        0.62.0   2022-04-19 [1] CRAN (R 4.1.2)
    ##  memoise            2.0.1    2021-11-26 [1] CRAN (R 4.1.2)
    ##  mime               0.12     2021-09-28 [1] CRAN (R 4.1.2)
    ##  modelr             0.1.8    2020-05-19 [1] CRAN (R 4.1.2)
    ##  munsell            0.5.0    2018-06-12 [1] CRAN (R 4.1.2)
    ##  nnet               7.3-16   2021-05-03 [4] CRAN (R 4.0.5)
    ##  pillar             1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
    ##  pkgconfig          2.0.3    2019-09-22 [1] CRAN (R 4.1.2)
    ##  plyr               1.8.7    2022-03-24 [1] CRAN (R 4.1.2)
    ##  png                0.1-7    2013-12-03 [1] CRAN (R 4.1.2)
    ##  polyclip           1.10-0   2019-03-14 [1] CRAN (R 4.1.2)
    ##  promises           1.2.0.1  2021-02-11 [1] CRAN (R 4.1.2)
    ##  purrr            * 0.3.4    2020-04-17 [1] CRAN (R 4.1.2)
    ##  R6                 2.5.1    2021-08-19 [1] CRAN (R 4.1.2)
    ##  RColorBrewer       1.1-3    2022-04-03 [1] CRAN (R 4.1.2)
    ##  Rcpp               1.0.9    2022-07-08 [1] CRAN (R 4.1.2)
    ##  RCurl              1.98-1.7 2022-06-09 [1] CRAN (R 4.1.2)
    ##  readr            * 2.1.2    2022-01-30 [1] CRAN (R 4.1.2)
    ##  readxl             1.4.0    2022-03-28 [1] CRAN (R 4.1.2)
    ##  reprex             2.0.1    2021-08-05 [1] CRAN (R 4.1.2)
    ##  rlang              1.0.4    2022-07-12 [1] CRAN (R 4.1.2)
    ##  rmarkdown          2.14     2022-04-25 [1] CRAN (R 4.1.2)
    ##  rpart              4.1-15   2019-04-12 [4] CRAN (R 4.0.0)
    ##  RSQLite            2.2.9    2021-12-06 [1] CRAN (R 4.1.2)
    ##  rstudioapi         0.13     2020-11-12 [1] CRAN (R 4.1.2)
    ##  rvest              1.0.2    2021-10-16 [1] CRAN (R 4.1.2)
    ##  S4Vectors          0.32.4   2022-03-24 [1] Bioconductor
    ##  scales             1.2.0    2022-04-13 [1] CRAN (R 4.1.2)
    ##  sessioninfo        1.2.2    2021-12-06 [1] CRAN (R 4.1.2)
    ##  shiny              1.7.1    2021-10-02 [1] CRAN (R 4.1.2)
    ##  stringi            1.7.6    2021-11-29 [1] CRAN (R 4.1.2)
    ##  stringr          * 1.4.0    2019-02-10 [1] CRAN (R 4.1.2)
    ##  survival           3.2-13   2021-08-24 [4] CRAN (R 4.1.1)
    ##  tibble           * 3.1.7    2022-05-03 [1] CRAN (R 4.1.2)
    ##  tidyr            * 1.2.0    2022-02-01 [1] CRAN (R 4.1.2)
    ##  tidyselect         1.1.2    2022-02-21 [1] CRAN (R 4.1.2)
    ##  tidyverse        * 1.3.1    2021-04-15 [1] CRAN (R 4.1.2)
    ##  tweenr             1.0.2    2021-03-23 [1] CRAN (R 4.1.2)
    ##  tzdb               0.3.0    2022-03-28 [1] CRAN (R 4.1.2)
    ##  UpSetR           * 1.4.0    2019-05-22 [1] CRAN (R 4.1.2)
    ##  utf8               1.2.2    2021-07-24 [1] CRAN (R 4.1.2)
    ##  vctrs              0.4.1    2022-04-13 [1] CRAN (R 4.1.2)
    ##  withr              2.5.0    2022-03-03 [1] CRAN (R 4.1.2)
    ##  xfun               0.31     2022-05-10 [1] CRAN (R 4.1.2)
    ##  XML                3.99-0.8 2021-09-17 [1] CRAN (R 4.1.2)
    ##  xml2               1.3.3    2021-11-30 [1] CRAN (R 4.1.2)
    ##  xtable             1.8-4    2019-04-21 [1] CRAN (R 4.1.2)
    ##  XVector            0.34.0   2021-10-26 [1] Bioconductor
    ##  yaml               2.3.5    2022-02-21 [1] CRAN (R 4.1.2)
    ##  zlibbioc           1.40.0   2021-10-26 [1] Bioconductor
    ## 
    ##  [1] /home/thomaz/R/x86_64-pc-linux-gnu-library/4.1
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
