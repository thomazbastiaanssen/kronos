Please note this package is under active construction.

# Kronos

If you use this software, please cite our work:

*\#include citation here with DOI*

The following document is adapted from the supplementary materials for
this manuscript.

## 0. Introduction

Here, we will demonstrate how to use the Kronos package to assess
circadian rhythms in biological data collected over the day. For this
demonstration we have adapted some data from our laboratory, currently
undergoing peer review (reference to come).

We will demonstrate three common examples of experimental design
currently encountered in circadian rhythms analysis:

1.  Circadian rhythm analysis of a single variable over a 24-hour period

2.  Circadian rhythm analysis of a single variable over a 24-hour period
    with two or more treatment groups

3.  Circadian rhythm analysis of omics (higher-dimensional) data over a
    24-hour period with two or more treatment groups

All R code used to transform, reorganise and plot the data is also shown
below to provide a toolkit for aspiring and veteran bioinformaticians
alike. It should be noted that some variables in the demonstration data
sets provided here have been manipulated to better demonstrate the
functionality of this package.

*At the end of this document, we have included a few excursions into
more advanced subjects we find useful, but that did not necessarily fit
in with the mainline analysis.*

### Code chunk: Install Kronos and load our libraries

``` r
#Install kronos
library(devtools)
#install_github("thomazbastiaanssen/kronos") 

#Load relevant packages
library(kronos)
library(tidyverse)
library(ggplot2)

#Load prepared data stored in Kronos library
data("kronos_demo")
```

### Code chunk: Load our data tables

Data should be prepared in “long form” for use in Kronos; that is with
values repeating in the “Timepoint” column, which defines when data was
collected during the period, here a 24-hour cycle.

Our package includes example datasets that we will use in this tutorial
that are pre-formatted. You can rearrange your data into long form using
pivot_longer() or gather() from the tidyverse.

Since we’re using prepared data, we already loaded it using
`data(kronos_demo)`. You can see how an example of how the data is
prepared below:

``` r
head(groupdata)
```

    ##   Animal_ID Timepoint Treatment Variable_1
    ## 1         6         5         A  1.2789856
    ## 2         7         5         A  1.0606877
    ## 3         8         5         A  1.0497167
    ## 4         9         5         A  1.0533610
    ## 5        10         5         A  1.4590203
    ## 6        12         5         A  0.8408964

Here we have prepared the omics data set with a separate metadata file
as is common when working with omics data sets. A metadata file can be
generated using select() in tidyr, or metadata and data can be combined
using inner_join if they contain an identical column (both name and
contents). The omics data set used here has been central-log transformed
to account for its compositional nature (see our guide
<https://arxiv.org/abs/2207.12475> for easy centred-log transformation).

## 1. Analysing Rhythmicity in a Single Group

We will start with the most simple example: analysing circadian
rhythmicity in a single experimental group for one outcome variable of
interest. For this we use the kronos function:

``` r
output <- kronos(formula = Variable_1 ~ time(Timepoint), 
                 data = onevariable, period = 24, 
                 verbose = T, pairwise = F)
```

    ## [1] "Using the following model: Variable_1 ~ Timepoint_cos + Timepoint_sin"
    ## [1] "Using the following model: Variable_1 ~ (Timepoint_cos + Timepoint_sin)"

Here we use the formula `Outcome Variable ~ time(Time Variable)`, which
is the most simple model used by the kronos function. We specify the
period as 24 (this can be adjusted as appropriate for the data
analysed). By selecting `verbose=T`, you will be able to see the models
run by the kronos function: this becomes increasingly useful when you
run more complex models. Finally we select `pairwise=F` here, as there
are no groups to compare for differences in rhythms.

The kronos function returns a kronosOut object, containing several
pieces of data that can be accessed using handy ‘getter’ functions,
which we will describe below:

1). The *getKronos_input()* function fetches the data that the model is
based on, as well as the calculated cosine and sine components.

``` r
head(getKronos_input(output))
```

    ##   Variable_1 Timepoint unique_group Timepoint_cos Timepoint_sin
    ## 1 -0.4239923         5         TRUE      0.258819     0.9659258
    ## 2 -1.0311723         5         TRUE      0.258819     0.9659258
    ## 3 -0.8739002         5         TRUE      0.258819     0.9659258
    ## 4 -0.5896825         5         TRUE      0.258819     0.9659258
    ## 5 -0.4174538         5         TRUE      0.258819     0.9659258
    ## 6 -0.4052512         5         TRUE      0.258819     0.9659258

2). The *getKronos_fit()* function fetches the key details for the
generated model that may be useful for prediction, modelling and other
statistical applications.

``` r
getKronos_fit(output)
```

    ## 
    ## Call:
    ## lm(formula = formula, data = data)
    ## 
    ## Coefficients:
    ##   (Intercept)  Timepoint_cos  Timepoint_sin  
    ##       -0.4063         0.3914        -0.6886

3). The *getKronos_trace()* function returns all the data required for
graphing the sinusoid curve, which can either be used in our specialized
ggplot2 functions, or can be used in other graphing packages. The
*y_hat* column represents the predicted value of the outcome variable:
this is essential for plotting the predicted sinusoid curve.

``` r
head(getKronos_trace(output))
```

    ##   Timepoint Timepoint_cos Timepoint_sin       y_hat unique_group
    ## 1      0.00     1.0000000    0.00000000 -0.01493121         TRUE
    ## 2      0.25     0.9978589    0.06540313 -0.06080489         TRUE
    ## 3      0.50     0.9914449    0.13052619 -0.10815799         TRUE
    ## 4      0.75     0.9807853    0.19509032 -0.15678776         TRUE
    ## 5      1.00     0.9659258    0.25881905 -0.20648595         TRUE
    ## 6      1.25     0.9469301    0.32143947 -0.25703975         TRUE

4). The *getKronos_groupwise()* function arguably fetches the most
useful output: this provides us with the p-value (p.val) and proportion
of the variance in the data explained (r.sq) when we fit our sinusoid
curve. Additionally we obtain the acrophase (acro) and amplitude of the
predicted curve, which can be used in our graphics functions to
visualise changes in curve with interventions (see `gg_kronos_circle`,
explored further below).

``` r
getKronos_groupwise(output)
```

    ##   unique_group        p.val      r.sq        avg     acro amplitude
    ## 1         TRUE 2.237134e-05 0.5345906 -0.4022588 19.97413  0.792033

### Figures

The package contains custom ggplot2 figure functions, that utilise the
kronos output to rapidly produce figures that convey important
information for circadian rhythms:

1.  `gg_kronos_circle()` generates a plot showing the acrophase and
    amplitude of the predicted curve, allowing the reader to rapidly
    access summary data regarding variables of interest, and to compare
    the summary data between groups in more complex models.At baseline,
    non-significant outcome measures are presented using dashed lines.
2.  `gg_kronos_sinusoid()` generates a x-y plot showing the outcome
    variable across the defined period. These graphs are useful for
    visualising the differences between specific timepoints assessed.

``` r
gg_kronos_circle(output)
```

![](README_files/figure-gfm/figures-1.png)<!-- -->

``` r
gg_kronos_sinusoid(output)
```

![](README_files/figure-gfm/figures-2.png)<!-- -->

## 2. Comparing Rhythmicity for More than Two Groups

Next we will demonstrate one of the unique features of the Kronos
package: the ability to compare circadian rhythms between more than two
groups. This is increasingly important as the use of complex
experimental designs grows in biological science. This example comprises
of three independent groups and is similar in setup to a one-way ANOVA.
For examples of more complex designs, see Excursion 1.

``` r
output2 <- kronos(formula = Variable_1 ~ Treatment + time(Timepoint), 
                  data = groupdata, period = 24, 
                  verbose = T, pairwise = T)
```

    ## [1] "Using the following model: Variable_1 ~ Treatment + Timepoint_cos + Timepoint_sin + Treatment:Timepoint_cos +     Treatment:Timepoint_sin - 1"
    ## [1] "Using the following model: Variable_1 ~ (Timepoint_cos + Timepoint_sin)"
    ## [1] "Using the following model: Variable_1 ~ (Timepoint_cos + Timepoint_sin)"
    ## [1] "Using the following model: Variable_1 ~ (Timepoint_cos + Timepoint_sin)"
    ## [1] "Fitting pairwise models"
    ## [1] "Using the following model: Variable_1 ~ unique_group * (Timepoint_cos + Timepoint_sin)"
    ## [1] "Using the following model: Variable_1 ~ unique_group * (Timepoint_cos + Timepoint_sin)"
    ## [1] "Using the following model: Variable_1 ~ unique_group * (Timepoint_cos + Timepoint_sin)"

``` r
gg_kronos_circle(output2)
```

![](README_files/figure-gfm/figures_complex-1.png)<!-- -->

``` r
gg_kronos_sinusoid(output2)
```

![](README_files/figure-gfm/figures_complex-2.png)<!-- -->

There are a few changes to the output generated by the kronos function:

1.  The *getKronos_groupwise()* output now contains a line for each of
    our groups.

``` r
getKronos_groupwise(output2)
```

    ##   unique_group       p.val       r.sq      avg      acro  amplitude
    ## 1            A 0.031490079 0.21886392 1.099699 12.681641 0.17966767
    ## 2            B 0.002912094 0.33147505 1.246347 18.592663 0.30283851
    ## 3            C 0.723230482 0.02287902 1.337834  3.447684 0.05150576

In this example you can see that groups A and B exhibit statistically
significant rhythms, while the model fitted to group C is
non-significant.

2.  We can now generate *pairwise models* output using `pairwise=T`.
    This generates pairwise comparisons between each of the groups:

``` r
getKronos_pairwise(output2)
```

    ## $`A vs B`
    ## Analysis of Variance Table
    ## 
    ## Response: Variable_1
    ##                            Df Sum Sq Mean Sq F value  Pr(>F)   
    ## unique_group                1 0.3626 0.36262  4.3279 0.04200 * 
    ## Timepoint_cos               1 0.1323 0.13228  1.5788 0.21406   
    ## Timepoint_sin               1 0.8859 0.88585 10.5728 0.00193 **
    ## unique_group:Timepoint_cos  1 0.4090 0.40901  4.8816 0.03118 * 
    ## unique_group:Timepoint_sin  1 0.5492 0.54916  6.5543 0.01314 * 
    ## Residuals                  57 4.7758 0.08379                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`A vs C`
    ## Analysis of Variance Table
    ## 
    ## Response: Variable_1
    ##                            Df Sum Sq Mean Sq F value    Pr(>F)    
    ## unique_group                1 0.9067 0.90674 14.3827 0.0003677 ***
    ## Timepoint_cos               1 0.1669 0.16692  2.6476 0.1093177    
    ## Timepoint_sin               1 0.0007 0.00068  0.0108 0.9174232    
    ## unique_group:Timepoint_cos  1 0.3425 0.34246  5.4321 0.0233965 *  
    ## unique_group:Timepoint_sin  1 0.0390 0.03899  0.6185 0.4349195    
    ## Residuals                  56 3.5305 0.06304                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`B vs C`
    ## Analysis of Variance Table
    ## 
    ## Response: Variable_1
    ##                            Df Sum Sq Mean Sq F value  Pr(>F)   
    ## unique_group                1 0.1279 0.12786  1.5595 0.21685   
    ## Timepoint_cos               1 0.0487 0.04868  0.5937 0.44417   
    ## Timepoint_sin               1 0.5623 0.56230  6.8581 0.01129 * 
    ## unique_group:Timepoint_cos  1 0.0025 0.00251  0.0306 0.86179   
    ## unique_group:Timepoint_sin  1 0.8940 0.89402 10.9038 0.00166 **
    ## Residuals                  57 4.6735 0.08199                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Above we can see that overall group A is significantly different between
B and C, and that group B exhibits a significantly different rhythm from
A and C.

3.  When including independent variables, kronos will also calculate an
    overall interaction with time which can be accessed using the code
    below:

``` r
getKronos_pairwise_p(output2)
```

    ##             adj.p.val
    ## Treatment 0.004532107

This is calculated by performing a Bonferroni correction on the
interactions between both the sine and cosine time components and the
independent variable. The p-value reported is the lowest following
correction.

## 3. ’Omics Analysis

Now we will demonstrate how to adapt the kronos package for ’omics
analysis, where there are many outcome variables. We have written the
`fw_kronos` function specifically for this purpose. This function
behaves very similar to the main kronos function.

It requires two core types of data: First, a table of data with rows as
features and columns as samples as input. Make sure that the feature
labels are row names rather than a column. Second, a metadata table with
rows as samples and metadata entries as columns. Suitable input data
looks like this:

``` r
head(bigdata)
```

    ##                   X49        X51        X52        X24        X25        X26
    ## Variable_1  0.6997887  0.7510796  0.4496661  0.7070641  0.7314069  0.4365048
    ## Variable_2  0.9185654  0.3832566  0.5772013  0.2469332  1.1264888  1.6429999
    ## Variable_3  5.7376737  5.7759893  5.8397621  4.6582202  4.6705684  4.7181954
    ## Variable_4 -0.3160217  0.2157887  5.9146992  0.1599933  0.2605589  7.5024676
    ## Variable_5 -4.8665730 -4.5445920 -3.5170020 -4.6925208 -4.3942593  3.3433245
    ## Variable_6 -0.9003465 -0.9120022 -1.0986671 -0.8211809 -1.3099717 -0.9428444
    ##                   X27        X28        X32        X36        X41        X42
    ## Variable_1  0.3050582  0.5632826  0.7671778  0.2908497  0.4404065  0.6242118
    ## Variable_2  0.8520508  0.5217500  0.8284091  1.1657981  1.1946169  0.8302731
    ## Variable_3  5.4136690  5.6849934  3.3614320  4.9325805  3.8249809  4.3401531
    ## Variable_4  7.3581878  0.2856433  7.7534970  7.3373731  6.2550336 -0.8893809
    ## Variable_5  3.3497446 -4.6677011 -4.9639562  2.7522326  2.3746696 -5.1336354
    ## Variable_6 -1.1136502 -1.0716995 -1.0774725 -1.4150184 -1.1520411 -1.0260923
    ##                   X43        X45        X47         X48        X91        X92
    ## Variable_1  0.7110339  0.8907024  0.5376684  0.65129006  0.6168283  0.5400244
    ## Variable_2  0.3342391  0.3322208  0.5134482  0.49903770 -0.2036331  1.1558801
    ## Variable_3  3.5427761  4.9931600  5.0252634  4.08640012  4.6361764  3.8853670
    ## Variable_4 -0.2464124  0.2829712 -0.5588745  0.02922933  0.1884404  6.9909271
    ## Variable_5 -5.0603399 -4.9920546 -4.8408525 -4.49220594 -4.8216887  2.7678792
    ## Variable_6 -1.0498812 -0.9101319 -1.1257748 -0.97311587 -0.9292009 -1.6540957
    ##                   X93        X94         X95        X96        X97        X98
    ## Variable_1  0.3100736  1.0005642  0.72362110  0.4967626  0.2618956  0.8005609
    ## Variable_2  0.5399346  0.6177746  0.72405229  1.4073206  1.1602874  0.5113400
    ## Variable_3  3.7830854  4.7669080  4.18681741  3.3145823  3.2333298  4.9457181
    ## Variable_4  7.0053741  0.2228974  0.02462879  7.6919773  7.5951861 -0.1985578
    ## Variable_5  2.1496092 -4.3879434 -4.63576141  2.2897306  2.4434397 -5.1091632
    ## Variable_6 -1.5787430 -1.4486857 -1.60453009 -1.2501476 -1.7315140 -0.9676350
    ##                   X54        X55         X56        X57        X59        X60
    ## Variable_1  0.5918080  0.7624666  0.74505297  0.6384300  0.5873153  0.3148437
    ## Variable_2  0.2558126  0.7986753  0.88080266  0.7036921  0.8342686  0.5704253
    ## Variable_3  5.8389649  6.0499947  6.13968090  6.2819807  5.4781121  6.2958542
    ## Variable_4 -0.3687695 -0.7125948 -0.04492208 -0.2076046  7.0065307  6.5890896
    ## Variable_5 -5.0480254 -5.4916166 -4.82088337 -4.6283552 -2.8305170  2.6337269
    ## Variable_6 -0.9324417 -0.7478880 -0.80980389 -0.8376803 -1.2139439 -0.8189141
    ##                   X61        X53         X6         X7         X8         X9
    ## Variable_1  0.3936574  0.6690287  0.4582703  0.3103756  0.3522794  0.7901186
    ## Variable_2  0.6329971  0.5440797  0.4955525  1.3053250  0.9315610  0.5548109
    ## Variable_3  6.3613095  6.2994679  5.0091207  5.3518483  4.4639989  5.5650364
    ## Variable_4  7.0880150 -0.1761231  6.4618090  7.1714075  6.8855375 -0.5202781
    ## Variable_5  2.9440566 -4.9637119 -4.2328681  3.5463528  3.0479420 -5.2488511
    ## Variable_6 -0.8510420 -0.7938576 -1.1617073 -1.1856011 -1.3188788 -0.9779349
    ##                   X10        X11        X12        X13        X16        X17
    ## Variable_1  0.7957730  0.9008205  0.6991660  0.7147391  0.7205021  0.2611259
    ## Variable_2  0.7497420  1.0883964  0.5020052  0.8623026  0.6734388  1.0862295
    ## Variable_3  4.5781266  4.0098843  4.6558234  5.3207146  5.1625729  3.8779918
    ## Variable_4 -0.2709393 -0.3357453 -0.4618575 -0.4217701 -0.1509541  7.1560720
    ## Variable_5 -5.1850036 -5.5402693 -4.7310010 -5.3859795 -5.1736359  3.0596881
    ## Variable_6 -1.2913439 -0.9627335 -0.8969915 -1.0836341 -1.1946675 -1.7234933
    ##                   X18        X19        X20       X101       X102        X66
    ## Variable_1  0.2293843  0.7105868  0.2259793  0.1652803  0.2696269  0.6140433
    ## Variable_2  0.8557440  0.5047304  1.2554064  1.3156165  1.1508952  0.3800611
    ## Variable_3  4.9297066  4.8625317  6.1120729  4.9211852  5.4037495  6.1319458
    ## Variable_4  7.3223489 -0.2968163  0.8280426  6.7753087  7.0185169 -0.1493441
    ## Variable_5  2.6618657 -5.2913340 -5.0843968  2.5483896  2.4561961 -6.8722489
    ## Variable_6 -1.3825037 -1.1444122 -1.4095839 -1.4145341 -0.9718657 -1.5317902
    ##                   X67        X68        X69        X70        X71       X122
    ## Variable_1  0.7456147  0.5829598  0.6087710  0.6079646  0.7151732  0.6890556
    ## Variable_2  0.2593484  0.3630193  0.4006584  0.4442866  0.5630537  0.9222034
    ## Variable_3  5.3555927  5.4054331  4.4615611  4.8859011  4.5141747  3.9319797
    ## Variable_4 -0.1575398 -0.4175065 -0.5709533  0.1754378 -0.3499401  0.4108061
    ## Variable_5 -4.1255322 -4.6299045 -4.8427934 -4.4781669 -4.9205045 -6.3323443
    ## Variable_6 -1.4938874 -1.2668969 -1.4396786 -1.2833023 -1.3753133 -1.0519939
    ##                  X123        X78        X79        X80        X81        X82
    ## Variable_1  0.7557379  0.5140498  0.7832171  0.7020577  0.7580827  0.7347557
    ## Variable_2  0.6811475  0.4506265  0.3478668  0.5012985  0.5267755  0.4573858
    ## Variable_3  3.9071926  4.9206138  4.6098578  4.4604895  4.1818153  3.9260191
    ## Variable_4  0.4005497 -0.1681288 -0.3665230 -0.7452182 -0.3277716  0.2835825
    ## Variable_5 -5.8645424 -4.4419256 -4.6515637 -6.8825993 -6.4151672 -4.2648530
    ## Variable_6 -1.1582234 -1.2806064 -1.0971614 -1.0648656 -1.0601368 -1.1149639
    ##                   X83       X117       X118         X72          X73
    ## Variable_1  0.4667721  0.5233439  0.6981221  0.28917763  0.545883127
    ## Variable_2  0.1116239  1.2661391  0.3044250  0.45741658  0.225358942
    ## Variable_3  4.7520230  5.8267016  5.6011169  4.43326051  3.868266095
    ## Variable_4  0.5045038  2.3221394  0.1970564  0.04569196  0.000308751
    ## Variable_5 -4.1757361 -4.9286835 -4.3020809 -4.25257028 -4.191003117
    ## Variable_6 -1.1149869 -1.8674129 -1.1634476 -1.22139052 -1.312496130
    ##                    X74        X75        X76         X77        X84        X86
    ## Variable_1  0.71646946  0.6637308  0.6331902  0.67592799  0.5891517  0.5634181
    ## Variable_2  0.36533193  0.5220832  0.4850906  0.38076857  0.6551901  0.6027763
    ## Variable_3  4.43997097  3.4729502  3.6269892  3.11411938  5.0924947  4.0079374
    ## Variable_4 -0.05991425 -0.3742904  5.1834396 -0.01419927  0.1253499  0.6806781
    ## Variable_5 -5.84929241 -5.8952158 -4.0259951 -4.35832164 -5.6764097 -3.2125444
    ## Variable_6 -1.25589978 -1.5256563 -1.4646955 -1.37751046 -1.1983062 -1.1447081
    ##                   X62        X63        X64        X65        X87        X88
    ## Variable_1  0.6169733  0.5554572  0.5879815  0.6416063  0.6990261  0.5855905
    ## Variable_2  0.3746583  0.6446069  0.4081766  0.7163491  0.4934286  1.1794391
    ## Variable_3  4.5038658  3.9332601  3.8058852  4.2996411  3.9030516  4.0056343
    ## Variable_4 -0.2730551 -0.2193646 -0.2588940 -0.1600187 -0.1442005  1.0764241
    ## Variable_5 -4.5504286 -5.1126163 -4.9608479 -5.4398616 -4.9390099 -6.0991683
    ## Variable_6 -1.0691214 -1.1858376 -0.9931351 -1.2324469 -0.9611707 -1.3540819
    ##                   X89        X90       X104       X105       X106       X107
    ## Variable_1  0.2557189  0.5945848  0.7972229  0.7564802  0.4739431  0.6940986
    ## Variable_2  0.9344481  0.4940206  0.8430432  0.3743948  0.7450903  0.9318353
    ## Variable_3  2.2381372  2.8617828  4.3390711  3.5502968  2.7524139  3.2552007
    ## Variable_4  1.4299143  1.2137917  0.3843475 -0.3697616  3.4090403 -0.5189534
    ## Variable_5  1.8578313 -1.3871499 -6.1177815 -5.6651668  1.5531021 -6.4803840
    ## Variable_6 -1.4206713 -1.1817701 -1.0585206 -1.2068126 -1.4310450 -1.0401184
    ##                  X108       X119       X121       X109       X110       X111
    ## Variable_1  0.6709748  0.4852421  0.5560595  0.5617770  0.5037281  0.5826479
    ## Variable_2  0.6184588  0.6118202  0.4687804  0.5137546  0.3432362  0.3546838
    ## Variable_3  3.4680358  4.5948820  3.7016101  3.7894121  3.3425325  3.8948287
    ## Variable_4 -0.3402560 -0.2216749  7.1128614 -0.0288680 -0.5828558 -0.4048614
    ## Variable_5 -5.2596745 -4.9723846 -4.6664862 -5.0187888 -5.1064492 -4.9731973
    ## Variable_6 -1.0332158 -1.3102901 -1.4175412 -1.1785680 -1.2481954 -0.9879016
    ##                  X112       X113       X114       X115       X116
    ## Variable_1  0.5442134  0.6115724  0.6995537  0.6529442  0.3673309
    ## Variable_2  0.7881500  1.0268203  0.5297310  0.2029693  1.1657880
    ## Variable_3  2.9838677  3.2566461  3.7821337  3.8743885  3.0907792
    ## Variable_4  1.8846074  2.1434427  7.5134831 -0.5780600  1.9926192
    ## Variable_5 -0.1634759  0.1940087 -4.5849179 -4.8306453 -2.9008829
    ## Variable_6 -1.1988142 -1.1683957 -1.0639460 -1.2451126 -1.0458468

``` r
head(bigmeta)
```

    ##   Animal_ID Group Timepoint
    ## 1       X49     B         5
    ## 2       X51     B         5
    ## 3       X52     B         5
    ## 4       X24     B         5
    ## 5       X25     B         5
    ## 6       X26     B         5

``` r
out_list = fw_kronos(bigdata, 
                     formula = ~ Group + time(Timepoint), 
                     metadata = bigmeta, 
                     pairwise = T) 
```

Now we have a list of kronosOut objects, which contain all our results.
This can be cumbersome to do manually, so we wrote the
`kronosListToTable` function for this purpose:

``` r
fit_df = kronosListToTable(out_list)

write.csv(fit_df, "README_files/RhythmicityResults.csv")
```

This will generate a csv containing the individual rhythmicity
calculations for the whole data set with an FDR correction to account
for multiple tests.

The resulting csv can be found
[here](https://github.com/thomazbastiaanssen/kronos/tree/main/README_files/RhythmicityResults.csv).

``` r
head(fit_df)
```

    ##                 B_p.val     A_p.val   C_p.val      B_r.sq     A_r.sq
    ## Variable_1 9.532798e-01 0.030901063 0.7771550 0.003294345 0.21991674
    ## Variable_2 9.297825e-01 0.002684955 0.3122346 0.005008424 0.34483205
    ## Variable_3 2.784019e-07 0.755935357 0.4488660 0.646891288 0.01978728
    ## Variable_4 8.603783e-01 0.009861269 0.4818005 0.010317651 0.28103213
    ## Variable_5 7.528003e-01 0.012385930 0.9204081 0.019392619 0.26923016
    ## Variable_6 3.588405e-03 0.103475111 0.3808739 0.321777005 0.14958456
    ##                 C_r.sq      B_avg      A_avg      C_avg    B_acro    A_acro
    ## Variable_1 0.017847065  0.5939175  0.5809418  0.5868535 12.787074 22.474987
    ## Variable_2 0.079780343  0.7216686  0.6931953  0.5981656  8.550615  8.798428
    ## Variable_3 0.055610401  4.9514601  4.8769289  3.7125717  1.243837 10.316266
    ## Variable_4 0.050822008  2.8118776  1.6005672  0.9716583 16.630514  8.855465
    ## Variable_5 0.005906638 -2.3709104 -3.5201319 -4.0940532 17.758696  9.103868
    ## Variable_6 0.066625775 -1.0943657 -1.2422211 -1.2090199  3.004538  3.396861
    ##               C_acro B_amplitude A_amplitude C_amplitude      B_q.val
    ## Variable_1 16.862949  0.01475693   0.1358320  0.02294777 9.982317e-01
    ## Variable_2 16.245056  0.03722544   0.2839437  0.10185205 9.982317e-01
    ## Variable_3  6.883162  1.09917381   0.1266424  0.20292811 1.559051e-05
    ## Variable_4 21.061482  0.51111164   2.2500130  0.66328131 9.982317e-01
    ## Variable_5 21.000586  0.68124941   2.4062525  0.23430428 9.982317e-01
    ## Variable_6 18.143617  0.21593906   0.1178712  0.05540318 6.698355e-02
    ##               A_q.val   C_q.val Group:Time_p.val Group:Time_q.val
    ## Variable_1 0.11536397 0.9909656     8.892379e-02     5.241823e-01
    ## Variable_2 0.03651062 0.7306471     6.412060e-02     4.824477e-01
    ## Variable_3 0.77508921 0.7935538     1.141716e-06     6.393608e-05
    ## Variable_4 0.07520935 0.7935538     1.376360e-01     5.246527e-01
    ## Variable_5 0.07520935 0.9937105     1.494553e-01     5.246527e-01
    ## Variable_6 0.22729268 0.7668614     3.274827e-02     4.548194e-01

We can use a similar approach to obtain other components of the
kronosOut objects. Below we include the code to obtain the pairwise
comparisons as a single csv, as this is slightly more difficult.

``` r
#Create an empty container list of the appropriate length
pairwise_list = vector(mode = "list", length = length(out_list)) 

#The for-loop below generates a list containing the pairwise test results
for(m in 1:length(out_list)){
  pairwise_list[[m]] <- out_list[[m]]@pairwise_models
  names(pairwise_list)[m] <- names(out_list)[m]
  
}

#Generate a single bound list
bound_list <- lapply(X = pairwise_list, FUN = function(x){do.call(rbind, x)})

#collapse the bound list for each variable into a single dataframe 
pairwise_df <- do.call(rbind, bound_list)

#separate the comparison and the effects to make results more readable
pairwise_csv = pairwise_df %>% 
  rownames_to_column("ID") %>%
  separate(col = ID, into = c("Feature","ID"), sep = "\\.", extra = "merge") %>%
  separate(col = ID, into = c("Comparison","Effect"), sep = "\\.")
  

write.csv(pairwise_csv, "README_files/PairwiseResults.csv") 
```

The resulting csv can be found
[here](https://github.com/thomazbastiaanssen/kronos/tree/main/README_files/PairwiseResults.csv).

### Figures

`gg_kronos_acrogram()` is a visualisation function we have designed
specifically for omics datasets. This function provides a polar
histogram of your dataset’s acrophases. This allows you to compare
overall rhythmicity between groups. Below we can see that a large
proportion of the variables peak between ZT20-23 in Group A, while the
variables in Groups B and C are less synchronous.

``` r
gg_kronos_acrogram(out_list)
```

![](README_files/figure-gfm/acrogram-1.png)<!-- -->

We can also use automation to obtain individual plots for our omics data
set. Here we will demonstrate how to obtain sinusoid curves for each
outcome measure in the data set.

``` r
 #Create an empty container list of the appropriate length
plot_list <- vector(mode = "list", length = length(out_list))

for(q in 1:length(out_list)){
  
  #save plot into relevant position in list
  plot_list[[q]] <- gg_kronos_sinusoid(out_list[[q]]) 
  
}

#to plot & save the feature graphs to a pdf:
pdf("README_files/plots_circadian.pdf")
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
invisible(dev.off())
```

The resulting pdf can be found
[here](https://github.com/thomazbastiaanssen/kronos/tree/main/README_files/plots_circadian.pdf).
The same approach can be applied for obtaining individual circle
figures.

## 4. Discussion

Here we have presented standard data for the analysis of time-of-day.
Some points to consider for your data is whether you can assume a
24-hour period, and whether your data is evenly distributed over the
period. Please note that you will require a minimum of three data points
over your period to make use of these functions, and indeed any function
using sinusoid models. Currently the kronos package is not able to
estimate period; a wide range of packages are capable of determining
your period if this is necessary for your research question. Please note
that period estimation requires even more temporal resolution: some
recommend a minimum of sampling every 2 hours over a 48-hr window
([Hughes et al.,
2017](https://journals.sagepub.com/doi/full/10.1177/0748730417728663)).

This tutorial is merely a template. Depending on your experimental
set-up, findings and experimental questions you may need to adjust your
approach. However, as complex statistical models become more frequent in
the study of circadian rhythms, functions that can incorporate more
complex design than two-group comparisons are essential for advancement
of the field.

We have provided figure generation functions as clear communication of
results is essential to producing good and useful science. Please see
below for more details for figure customisation. We hope that both
aspiring and veteran bioinformaticians and circadian rhythm biologists
will find our guide helpful.

## Excursion 1. Customising Figures

The two figure functions, `gg_kronos_circle()` and
`gg_kronos_sinusoid()`, are designed to be fully compatible with ggplot2
and therefore are fully customisable. Below is an example of how one
could customise kronos plots using ggplot2 syntax

``` r
gg_kronos_circle(output2) +
  scale_fill_manual(values = c("A" = "#169B62",
                               "B" = "#FFFFFF", 
                               "C" = "#FF883E")) +
  ggtitle("Figure title")
```

![](README_files/figure-gfm/customise%20figure-1.png)<!-- -->

We encourage users to take advantage of the extensive range of online
tutorials and add-on packages available for ggplot2.

## Excursion 2. Assessing Rhythmicity in More Complex Models

One of the key features of this package is the use of a formula input,
which allows for analysis of complex models. Below we will demonstrate
how kronos performs when assessing data with two independent variables.

``` r
data3 <- twowaydata

two.way.data.long <- data3 %>%
  pivot_longer(cols=starts_with("Variable_"), names_to = "Variables", values_to = "Value") %>%
  as.data.frame()

#collect all outcome variables
two_way_data_names <- unique(two.way.data.long$Variables) 

#Create an empty container list of the appropriate length
data.list <- vector(mode = "list", length = length(two_way_data_names)) 

for(n in 1:length(two_way_data_names)){
  data.list[[n]] <- two.way.data.long %>% filter(Variables == two_way_data_names[n])
}

two_way_out_list = lapply(X = data.list, 
                          FUN = function(y){kronos(data = y, 
                                                   Value ~ Factor_A*Factor_B + time(Timepoint), 
                                                   period = 24, pairwise = T, verbose = F)
                            }
                          )

names(two_way_out_list) <- two_way_data_names


gg_kronos_sinusoid(two_way_out_list$Variable_1)
```

![](README_files/figure-gfm/two-way%20output-1.png)<!-- -->

``` r
gg_kronos_sinusoid(two_way_out_list$Variable_2)
```

![](README_files/figure-gfm/two-way%20output-2.png)<!-- -->

``` r
gg_kronos_sinusoid(two_way_out_list$Variable_3)
```

![](README_files/figure-gfm/two-way%20output-3.png)<!-- -->

``` r
gg_kronos_sinusoid(two_way_out_list$Variable_4)
```

![](README_files/figure-gfm/two-way%20output-4.png)<!-- -->

Here we have analysed 4 outcome variables which all show different
interaction effects. Here we will go into depth examining the effects
observed in Variable 1 as an example of how to interpret kronos output
for more complex designs.

``` r
gg_kronos_sinusoid(two_way_out_list$Variable_1)
```

![](README_files/figure-gfm/two-way%20var1%20graphs-1.png)<!-- -->

``` r
gg_kronos_circle(two_way_out_list$Variable_1)
```

![](README_files/figure-gfm/two-way%20var1%20graphs-2.png)<!-- -->

1). As before, we can use the `getKronos_groupwise()` function to obtain
individual rhythmicity for each group. Here you can see that 3/4 groups
exhibit rhythmicity and that both conventional groups share a similar
acrophase (which is illustrated in the figures above as well).

``` r
getKronos_groupwise(two_way_out_list$Variable_1)
```

    ##           unique_group       p.val      r.sq      avg      acro amplitude
    ## 1   Antibiotics_Stress 0.003290904 0.2855699 4712.202 16.051311 548.95960
    ## 2  Antibiotics_Control 0.045674193 0.1575633 2340.530  3.493713 219.93545
    ## 3  Conventional_Stress 0.001783400 0.3267086 2329.867 18.079296 256.07917
    ## 4 Conventional_Control 0.075062913 0.1452404 1349.590 18.002807  91.75282

2). With the `getKronos_pairwise_p()` function we can assess the
interaction of each of our experimental factors with the time component
of the model: here you can see that both main effects and the
interaction significantly interact with the time component.

``` r
getKronos_pairwise_p(two_way_out_list$Variable_1)
```

    ##                      adj.p.val
    ## Factor_A          1.000000e+00
    ## Factor_B          3.959353e-05
    ## Factor_A:Factor_B 1.911130e-02

3). Next we can use the `getKronos_pairwise()` function to obtain the
pairwise group comparisons. This allows us to determine how each group
differs from one another. For example, here you can see that
Conventional+Stress and Antibiotics+Control only exhibit a significant
group\*Timepoint_sin interaction. This is unsurprising as the groups
have the same average value but exhibit a rhythm shifted by 12 hours.

``` r
getKronos_pairwise(two_way_out_list$Variable_1)
```

    ## $`Antibiotics_Stress vs Antibiotics_Control`
    ## Analysis of Variance Table
    ## 
    ## Response: Value
    ##                            Df    Sum Sq   Mean Sq  F value    Pr(>F)    
    ## unique_group                1 107088781 107088781 392.5199 < 2.2e-16 ***
    ## Timepoint_cos               1    174790    174790   0.6407 0.4261766    
    ## Timepoint_sin               1    486557    486557   1.7834 0.1860572    
    ## unique_group:Timepoint_cos  1   1596022   1596022   5.8500 0.0181806 *  
    ## unique_group:Timepoint_sin  1   4198056   4198056  15.3874 0.0002021 ***
    ## Residuals                  70  19097668    272824                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Antibiotics_Stress vs Conventional_Stress`
    ## Analysis of Variance Table
    ## 
    ## Response: Value
    ##                            Df    Sum Sq   Mean Sq  F value    Pr(>F)    
    ## unique_group                1 102318455 102318455 420.7923 < 2.2e-16 ***
    ## Timepoint_cos               1    741756    741756   3.0505   0.08536 .  
    ## Timepoint_sin               1   4761304   4761304  19.5812 3.697e-05 ***
    ## unique_group:Timepoint_cos  1    672717    672717   2.7666   0.10099    
    ## unique_group:Timepoint_sin  1    452020    452020   1.8590   0.17738    
    ## Residuals                  66  16048341    243157                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Antibiotics_Stress vs Conventional_Control`
    ## Analysis of Variance Table
    ## 
    ## Response: Value
    ##                            Df    Sum Sq   Mean Sq  F value    Pr(>F)    
    ## unique_group                1 210361830 210361830 969.9952 < 2.2e-16 ***
    ## Timepoint_cos               1    726171    726171   3.3484 0.0717181 .  
    ## Timepoint_sin               1   2785092   2785092  12.8423 0.0006375 ***
    ## unique_group:Timepoint_cos  1    669456    669456   3.0869 0.0834930 .  
    ## unique_group:Timepoint_sin  1   1403802   1403802   6.4730 0.0132640 *  
    ## Residuals                  67  14530220    216869                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Antibiotics_Control vs Conventional_Stress`
    ## Analysis of Variance Table
    ## 
    ## Response: Value
    ##                            Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## unique_group                1    1922    1922  0.0163 0.8988585    
    ## Timepoint_cos               1  168888  168888  1.4301 0.2358990    
    ## Timepoint_sin               1   18570   18570  0.1572 0.6929452    
    ## unique_group:Timepoint_cos  1  209481  209481  1.7739 0.1873534    
    ## unique_group:Timepoint_sin  1 1847703 1847703 15.6460 0.0001847 ***
    ## Residuals                  68 8030383  118094                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Antibiotics_Control vs Conventional_Control`
    ## Analysis of Variance Table
    ## 
    ## Response: Value
    ##                            Df   Sum Sq  Mean Sq  F value    Pr(>F)    
    ## unique_group                1 19495877 19495877 206.5666 < 2.2e-16 ***
    ## Timepoint_cos               1   175177   175177   1.8561  0.177510    
    ## Timepoint_sin               1   101312   101312   1.0734  0.303786    
    ## unique_group:Timepoint_cos  1   195432   195432   2.0707  0.154676    
    ## unique_group:Timepoint_sin  1   729445   729445   7.7288  0.006999 ** 
    ## Residuals                  69  6512262    94381                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Conventional_Stress vs Conventional_Control`
    ## Analysis of Variance Table
    ## 
    ## Response: Value
    ##                            Df   Sum Sq  Mean Sq  F value    Pr(>F)    
    ## unique_group                1 18114110 18114110 340.0055 < 2.2e-16 ***
    ## Timepoint_cos               1     1535     1535   0.0288   0.86575    
    ## Timepoint_sin               1  1116016  1116016  20.9479 2.184e-05 ***
    ## unique_group:Timepoint_cos  1      132      132   0.0025   0.96050    
    ## unique_group:Timepoint_sin  1   256055   256055   4.8062   0.03194 *  
    ## Residuals                  65  3462935    53276                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Session Info

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_IE.UTF-8        LC_COLLATE=en_GB.UTF-8    
    ##  [5] LC_MONETARY=en_IE.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    ##  [7] LC_PAPER=en_IE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4    
    ##  [5] readr_2.1.2     tidyr_1.2.0     tibble_3.1.7    ggplot2_3.3.6  
    ##  [9] tidyverse_1.3.1 kronos_0.1.0.0  devtools_2.4.3  usethis_2.1.5  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lubridate_1.8.0   prettyunits_1.1.1 ps_1.7.0          assertthat_0.2.1 
    ##  [5] rprojroot_2.0.3   digest_0.6.29     utf8_1.2.2        cellranger_1.1.0 
    ##  [9] R6_2.5.1          backports_1.4.1   reprex_2.0.1      evaluate_0.15    
    ## [13] highr_0.9         httr_1.4.3        pillar_1.7.0      rlang_1.0.4      
    ## [17] readxl_1.4.0      rstudioapi_0.13   callr_3.7.0       rmarkdown_2.14   
    ## [21] labeling_0.4.2    desc_1.4.1        munsell_0.5.0     broom_0.8.0      
    ## [25] compiler_4.1.2    modelr_0.1.8      xfun_0.31         pkgconfig_2.0.3  
    ## [29] pkgbuild_1.3.1    htmltools_0.5.2   tidyselect_1.1.2  fansi_1.0.3      
    ## [33] tzdb_0.3.0        crayon_1.5.1      dbplyr_2.2.0      withr_2.5.0      
    ## [37] brio_1.1.3        grid_4.1.2        jsonlite_1.8.0    gtable_0.3.0     
    ## [41] lifecycle_1.0.1   DBI_1.1.2         magrittr_2.0.3    scales_1.2.0     
    ## [45] cli_3.3.0         stringi_1.7.6     cachem_1.0.6      farver_2.1.1     
    ## [49] fs_1.5.2          remotes_2.4.2     testthat_3.1.4    xml2_1.3.3       
    ## [53] ellipsis_0.3.2    generics_0.1.2    vctrs_0.4.1       tools_4.1.2      
    ## [57] glue_1.6.2        hms_1.1.1         processx_3.6.0    pkgload_1.2.4    
    ## [61] fastmap_1.1.0     yaml_2.3.5        colorspace_2.0-3  sessioninfo_1.2.2
    ## [65] rvest_1.0.2       memoise_2.0.1     knitr_1.39        haven_2.5.0
