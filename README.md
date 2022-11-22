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

1.  Circadian rhythm analysis of 1 variable over a 24-hour period

2.  Circadian rhythm analysis of 1 variable over a 24-hour period with
    two or more treatment groups

3.  Circadian rhythm analysis of omics data over a 24-hour period with
    two or more treatment groups

All R code used to transform, reorganise and plot the data is also shown
below to hopefully provide a toolkit for aspiring and veteran
bioinformaticians alike. It should be noted that the analysis performed
here may not perfectly correspond to that performed in the published
manuscript: some variables in the datasets provided have been
manipulated to better demonstrate functions of the package.

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
###Thomaz to help wrangle this
```

### Code chunk: Load our data tables

Data should be prepared in “long form” for use in Kronos; that is with
values repeating in the “Timepoint” column, which defines when data was
collected in a 24-hour cycle.

Our package includes example datasets that we will use in this tutorial
that are pre-formatted. You can rearrange your data into long form using
pivot_longer() or gather() in tidyr.

Since we’re using prepared data, we already loaded it using
‘data*(INSERT DATA NAME HERE)*’ but typically we would do something like
this:

``` r
data("kronos_demo")
data1 <-   onevariable #One variable, one group
data2 <-   groupdata #One variable, more than two groups
bigdata <- bigdata #Omics dataset, more than two groups
bigmeta <- bigmeta #Omics sample metadata
```

You can see how an example of how the data is prepared below:

    ##   Animal_ID Timepoint Treatment Variable_1
    ## 1         6         5         A  1.2789856
    ## 2         7         5         A  1.0606877
    ## 3         8         5         A  1.0497167
    ## 4         9         5         A  1.0533610
    ## 5        10         5         A  1.4590203
    ## 6        12         5         A  0.8408964

Here we have prepared the omics dataset with a separate metadata file as
is common when working with omics datasets. A metadata file can be
generated using select() in tidyr, or metadata and data can be combined
using inner_join if they contain an identical column (both name and
contents). The omics dataset used here has been central-log transformed
to account for its compositional nature (see our guide
<https://arxiv.org/abs/2207.12475> for easy centred-log transformation).

## 1. Statistical Approach

*THOMAZ PLEASE WRITE SOMETHING*

## 2. Analysing Rhythmicity in a Single Group

We will start with the most simple example: analysing circadian
rhythmicity in a single experimental group for one outcome variable of
interest. For this we use the kronos function:

``` r
output <- kronos(formula = Variable_1 ~ time(Timepoint), 
                 data = data1, period = 24, 
                 verbose = T, pairwise = F)
```

    ## [1] "Using the following model: Variable_1 ~ Timepoint_cos + Timepoint_sin"
    ## [1] "Using the following model: Variable_1 ~ (Timepoint_cos + Timepoint_sin)"

Here we use the formula `Outcome Variable ~ time(Time Variable)`, which
is the most simple model used by the Kronos function. We specify the
period as 24 (this can be adjusted as appropriate for the data
analysed). By selecting `verbose=T`, you will be able to see the models
run by the kronos function: this becomes increasingly useful when you
run more complex models. Finally we select `pairwise=F` here, as there
are no groups to compare for differences in rhythms.

The kronos function creates a kronosOut object, containing several
lists, which we will describe below:

1). The *input* contains the data that the model is based on, as well as
the calculated cosine and sine components.

    ##   Variable_1 Timepoint unique_group Timepoint_cos Timepoint_sin
    ## 1 -0.4239923         5         TRUE      0.258819     0.9659258
    ## 2 -1.0311723         5         TRUE      0.258819     0.9659258
    ## 3 -0.8739002         5         TRUE      0.258819     0.9659258
    ## 4 -0.5896825         5         TRUE      0.258819     0.9659258
    ## 5 -0.4174538         5         TRUE      0.258819     0.9659258
    ## 6 -0.4052512         5         TRUE      0.258819     0.9659258

2). The *fit* contains the key coefficients for the model.

3). The *to_plot* contains all the data required for graphing the
sinusoid curve, which can either be used in our custom ggplot2
functions, or can be used in other graphing packages. *y_hat* is
predicted value of the outcome variable: this is essential for plotting
the predicted sinusoid curve.

    ##   Timepoint Timepoint_cos Timepoint_sin       y_hat unique_group
    ## 1      0.00     1.0000000    0.00000000 -0.01493121         TRUE
    ## 2      0.25     0.9978589    0.06540313 -0.06080489         TRUE
    ## 3      0.50     0.9914449    0.13052619 -0.10815799         TRUE
    ## 4      0.75     0.9807853    0.19509032 -0.15678776         TRUE
    ## 5      1.00     0.9659258    0.25881905 -0.20648595         TRUE
    ## 6      1.25     0.9469301    0.32143947 -0.25703975         TRUE

4). The *ind_fit* is arguably the most useful output: this provides us
with the p-value (p.val) and proportion of the variance in the data
explained (r.sq) when we fit our sinusoid curve. Additionally we obtain
the acrophase (acro) and amplitude of the predicted curve, which can be
used in our graphics functions to visualise changes in curve with
interventions (explored further in Sections 3 and 4).

    ##   unique_group        p.val      r.sq        avg     acro amplitude
    ## 1         TRUE 2.237134e-05 0.5345906 -0.4022588 19.97413  0.792033

### Figures

## 3. Comparing Rhythmicity for More than Two Groups

Next we will demonstrate one of the unique features of the Kronos
package: the ability to compare circadian rhythms between more than two
groups. This is increasingly important as the use of complex
experimental designs has increased in biological science. This example
comprises of three independent groups and is similar in setup to a
one-way ANOVA. For examples of more complex designs, see Excursion 1.

### Figures

## 4. Omics Analysis

BLEURGH

### Figures

## 5. Discussion

BLEURGH

## Excursion 1. Assessing Rhythmicity in More Complex Models
