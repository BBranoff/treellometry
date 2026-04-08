ecoinfromatics publication placeholder
================
Ben Branoff
April 8, 2026

# Coming Soon

This is the future site of the repository for the journal article: Price
et al. (2026) Local vs global allometric equations: Evaluating multiple
variables in mixed effects models for biomass estimation using a global
mangrove data compilation.

# Abstract

Accurate estimation of mangrove biomass is key component of the blue
carbon economy and allometric equations are a cornerstone of
methodological approaches from plot to landscape scales. A perennial
question facing investigators is whether global allometric equations can
be applied to their biome of interest, or whether local trees need to be
harvested to provide accurate estimates. Further, which combinations of
tree dimensions provide the most accurate estimates is not always clear.
We examined 16 different allometric equations for four tree measures
within a global data set for three prominent trans-Atlantic mangroves
species. We used mixed effects models with site location and tree
species modeled as random effects, and tree DBH, height, canopy diameter
and wood density as fixed effects. While the mixed effects model with
the greatest number of parameters was the best performing model based on
five different statistical measures with a mean error of 1.392 kg, the
improvement over simpler models was marginal. Simply measuring the DBH
of each tree would only add an additional 380 grams of error per tree
compared to the best performing model. These results provide qualified
support for the use of global equations and can help subsequent
investigators to prioritize data collection efforts.

# Methods

The below methods can be done manually, by downloading the referenced
data and functions from this repository and sourcing them in R, or they
can be done by installing the package in R.

``` r
##  to install this repository directly, the installation function (remotes) first need to be installed
install.packages("remotes")
##  then install from github
remotes::install_github("BBranoff/treellometry@ecoinf")
library(treellometry)
```

With the library successfully loaded, the following methods will
demonstrate the workflow documented in the above publication. The first
steps are to create the composite variables. These are commonly used in
allometric modelling because two allometric covariables are likely to be
correlated, violating the assumptions of multivariate models. Creating a
composite variable eliminates this issue, although it does introduce
interpretation considerations. Here, we also set the response variable,
which in our case is always aboveground biomass, as well as the various
predictor variables. We also establish the grouping variables, these
will be used to model random effects, which aims to investigate a
central premise of this publication, that species and location are
significant and substantial in influencing biomass predictions from
allometric measurements. In short, mixed effects models will ‘partition’
the variance in biomass associated with individual species and locations
and then attempt to model the remaining variance based on the allometric
predictor variables. This is very useful in determining how much of the
biomass is species or location specific, and how much depends primarily
on the allometric measurements. If a significant and substantial portion
of the variance in biomass can be attributed to these grouping
variables, that justifies the harvest of local trees in order to
determine location and species specific allometric models. On the other
hand, if the mixed effects models determine that an “insubstantial”
portion of the variance is due to these grouping variables, harvesting
local trees is probably not warranted and the use of ‘global’ models
that apply to all species and all locations is ‘good enough’.
Ultimately, that decision lies with local managers and scientists, but
here we provide the statistical information necessary to make that
informed decision. For a more detailed outline of mixed-effects
modelling in ecology, see [Harrison et
al. (2018)](https://doi.org/10.7717/peerj.4794).

## A look inside ‘explore_allom_models()’

The above use of the ‘explore_allom_models()’ accomplishes the bulk of
the analysis for this publication, which includes fitting models of
biomass from various combinations of predictor variables and fixed and
mixed effects, as well as retrieving model coefficients and fitness
metrics used in model evaluation. Rather than leaving all of that in the
black box of the function, here we break down the internal steps of
‘explore_allom_models()’ to provide more insight.

The first step is to check to make sure inputs are valid and then to run
two important transformations on the data. The first is a log
transformation, which is very common in allometry (see [Gingerich
(2000)](https://doi.org/10.1006/jtbi.2000.2008), and [Kerkhoff and
Enquist (2009)](https://doi.org/10.1016/j.jtbi.2008.12.026)). The second
is a scaling that is important for mixed-effects model performance and
interpretability [Harrison et
al. (2018)](https://doi.org/10.7717/peerj.4794). It is important to note
here, however, that transformations fundamentally change model
coefficients and their interpretation of the original data. For this
reason, coefficients for this publication are reported on models from
the non-scaled data.

``` r
  ##  check inputs
  if (responseVar %in% predictorVars) stop("Response variable found in predictor variables")
  # Log-transform response and predictors
  dat[[paste0("log", responseVar)]] <- log(dat[[responseVar]])
  for (var in predictorVars) {
    dat[[paste0("log", var)]] <- log(dat[[var]])
  }
  if (scle==TRUE){
    dat <- dat |>
      mutate( across(contains(predictorVars),function(x) scale(x)[,1]))#,.names = "{paste0(col, '_scaled')}"))
  }
  dat <- dat |> mutate( ID=1:n())
  predsMM <- predsF <- dat
  results <- list()
  coefs <- list()
  mods <- list()
  varGroup=model=0
```

Next, we begin to loop through each of the predictor variables and build
and fit models with different combinations if those variables. Here, we
use k to denote the number of predictor variables in the equation. A k
of 1 is a simple regression and a k greater than one is multiple
regression with more than one predictor. Two important exclusions are
happening below. The first is that the ‘WoodDensity’ variable is not
included in simple regression because its values are generalized for
each species and are not location specific. Second, is that composite
variables are not included in multiple regression because they are
already a combination of multiple variables and it simply doesnt make
sense in our case to include the influence of a variable twice in the
same model.

``` r
for (k in 1:length(predictorVars)) {
    if (k==1){
      ### create predictor combination sets
      ### exclude wood density because it is only from the San Juan dataset and only relative to entire species, not individuals
      combos <- combn(predictorVars[!grepl("WoodDensity",predictorVars)], k, simplify = FALSE)
    }else{
      ####  exclude the composite variables because they are already combined, doesnt make sense to include them in multivariate models
      if (grepl("Comp",predictorVars[k])){next}
      combos <- combn(predictorVars[!grepl("Comp",predictorVars)], k, simplify = FALSE)
    }
```

The ‘combos’ variable created in the above portion of the script is a
set of variable combinations that will be used in the model. They are
sent along to the next step, which further expands upon the ‘combos’ by
computing all possible combinations of those variables. This is
important because evaluating variable importance depends upon the order
in which a variable is present within a multiple regression model. As an
example, here is the result of combining all of the predictor variables
and then permutating their combinations

``` r
combos <- combn(predictorvars[!grepl("Comp",predictorvars)], 4, simplify = FALSE)
for (predss in combos) {
      ### establish all potential combinations of the current predictor set (change the order of predictors).
      perms <- gtools::permutations(length(predss),length(predss),predss)
      if (nrow(perms)==0){perms <- matrix(1,nrow = 1,ncol=1)}
}
print(perms)
```

    ##       [,1]                [,2]                [,3]               
    ##  [1,] "CanopyDiameter.m"  "DBH.cm"            "Height.m"         
    ##  [2,] "CanopyDiameter.m"  "DBH.cm"            "WoodDensity.g.cm3"
    ##  [3,] "CanopyDiameter.m"  "Height.m"          "DBH.cm"           
    ##  [4,] "CanopyDiameter.m"  "Height.m"          "WoodDensity.g.cm3"
    ##  [5,] "CanopyDiameter.m"  "WoodDensity.g.cm3" "DBH.cm"           
    ##  [6,] "CanopyDiameter.m"  "WoodDensity.g.cm3" "Height.m"         
    ##  [7,] "DBH.cm"            "CanopyDiameter.m"  "Height.m"         
    ##  [8,] "DBH.cm"            "CanopyDiameter.m"  "WoodDensity.g.cm3"
    ##  [9,] "DBH.cm"            "Height.m"          "CanopyDiameter.m" 
    ## [10,] "DBH.cm"            "Height.m"          "WoodDensity.g.cm3"
    ## [11,] "DBH.cm"            "WoodDensity.g.cm3" "CanopyDiameter.m" 
    ## [12,] "DBH.cm"            "WoodDensity.g.cm3" "Height.m"         
    ## [13,] "Height.m"          "CanopyDiameter.m"  "DBH.cm"           
    ## [14,] "Height.m"          "CanopyDiameter.m"  "WoodDensity.g.cm3"
    ## [15,] "Height.m"          "DBH.cm"            "CanopyDiameter.m" 
    ## [16,] "Height.m"          "DBH.cm"            "WoodDensity.g.cm3"
    ## [17,] "Height.m"          "WoodDensity.g.cm3" "CanopyDiameter.m" 
    ## [18,] "Height.m"          "WoodDensity.g.cm3" "DBH.cm"           
    ## [19,] "WoodDensity.g.cm3" "CanopyDiameter.m"  "DBH.cm"           
    ## [20,] "WoodDensity.g.cm3" "CanopyDiameter.m"  "Height.m"         
    ## [21,] "WoodDensity.g.cm3" "DBH.cm"            "CanopyDiameter.m" 
    ## [22,] "WoodDensity.g.cm3" "DBH.cm"            "Height.m"         
    ## [23,] "WoodDensity.g.cm3" "Height.m"          "CanopyDiameter.m" 
    ## [24,] "WoodDensity.g.cm3" "Height.m"          "DBH.cm"           
    ##       [,4]               
    ##  [1,] "WoodDensity.g.cm3"
    ##  [2,] "Height.m"         
    ##  [3,] "WoodDensity.g.cm3"
    ##  [4,] "DBH.cm"           
    ##  [5,] "Height.m"         
    ##  [6,] "DBH.cm"           
    ##  [7,] "WoodDensity.g.cm3"
    ##  [8,] "Height.m"         
    ##  [9,] "WoodDensity.g.cm3"
    ## [10,] "CanopyDiameter.m" 
    ## [11,] "Height.m"         
    ## [12,] "CanopyDiameter.m" 
    ## [13,] "WoodDensity.g.cm3"
    ## [14,] "DBH.cm"           
    ## [15,] "WoodDensity.g.cm3"
    ## [16,] "CanopyDiameter.m" 
    ## [17,] "DBH.cm"           
    ## [18,] "CanopyDiameter.m" 
    ## [19,] "Height.m"         
    ## [20,] "DBH.cm"           
    ## [21,] "Height.m"         
    ## [22,] "CanopyDiameter.m" 
    ## [23,] "DBH.cm"           
    ## [24,] "CanopyDiameter.m"

Next, each of the permutations in ‘perms’ is used in both a fixed
effects model as well as various mixed effects models.
