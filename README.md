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

``` r
data("mangroves")

###  create composite variables
mangroves <- mangroves |>
  mutate(Comp.DBH.H.Den = (DBH.cm^2)*Height.m*WoodDensity.g.cm3,
         Comp.Can.H.Den = (CanopyDiameter.m^2)*Height.m*WoodDensity.g.cm3)

###  set the variables of interest
responsevar <- "AGB.kg"
predictorvars <- c("DBH.cm","Height.m","CanopyDiameter.m","WoodDensity.g.cm3","Comp.DBH.H.Den","Comp.Can.H.Den")
###  the random effect variables
groupvars = c("Species", "Site")

mods <- explore_allom_models(mangroves,responsevar,predictorvars,groupvars)
```

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

    ##       [,1]                [,2]                [,3]                [,4]               
    ##  [1,] "CanopyDiameter.m"  "DBH.cm"            "Height.m"          "WoodDensity.g.cm3"
    ##  [2,] "CanopyDiameter.m"  "DBH.cm"            "WoodDensity.g.cm3" "Height.m"         
    ##  [3,] "CanopyDiameter.m"  "Height.m"          "DBH.cm"            "WoodDensity.g.cm3"
    ##  [4,] "CanopyDiameter.m"  "Height.m"          "WoodDensity.g.cm3" "DBH.cm"           
    ##  [5,] "CanopyDiameter.m"  "WoodDensity.g.cm3" "DBH.cm"            "Height.m"         
    ##  [6,] "CanopyDiameter.m"  "WoodDensity.g.cm3" "Height.m"          "DBH.cm"           
    ##  [7,] "DBH.cm"            "CanopyDiameter.m"  "Height.m"          "WoodDensity.g.cm3"
    ##  [8,] "DBH.cm"            "CanopyDiameter.m"  "WoodDensity.g.cm3" "Height.m"         
    ##  [9,] "DBH.cm"            "Height.m"          "CanopyDiameter.m"  "WoodDensity.g.cm3"
    ## [10,] "DBH.cm"            "Height.m"          "WoodDensity.g.cm3" "CanopyDiameter.m" 
    ## [11,] "DBH.cm"            "WoodDensity.g.cm3" "CanopyDiameter.m"  "Height.m"         
    ## [12,] "DBH.cm"            "WoodDensity.g.cm3" "Height.m"          "CanopyDiameter.m" 
    ## [13,] "Height.m"          "CanopyDiameter.m"  "DBH.cm"            "WoodDensity.g.cm3"
    ## [14,] "Height.m"          "CanopyDiameter.m"  "WoodDensity.g.cm3" "DBH.cm"           
    ## [15,] "Height.m"          "DBH.cm"            "CanopyDiameter.m"  "WoodDensity.g.cm3"
    ## [16,] "Height.m"          "DBH.cm"            "WoodDensity.g.cm3" "CanopyDiameter.m" 
    ## [17,] "Height.m"          "WoodDensity.g.cm3" "CanopyDiameter.m"  "DBH.cm"           
    ## [18,] "Height.m"          "WoodDensity.g.cm3" "DBH.cm"            "CanopyDiameter.m" 
    ## [19,] "WoodDensity.g.cm3" "CanopyDiameter.m"  "DBH.cm"            "Height.m"         
    ## [20,] "WoodDensity.g.cm3" "CanopyDiameter.m"  "Height.m"          "DBH.cm"           
    ## [21,] "WoodDensity.g.cm3" "DBH.cm"            "CanopyDiameter.m"  "Height.m"         
    ## [22,] "WoodDensity.g.cm3" "DBH.cm"            "Height.m"          "CanopyDiameter.m" 
    ## [23,] "WoodDensity.g.cm3" "Height.m"          "CanopyDiameter.m"  "DBH.cm"           
    ## [24,] "WoodDensity.g.cm3" "Height.m"          "DBH.cm"            "CanopyDiameter.m"

Next, each of the permutations in ‘perms’ (each row of the above table)
is used in both a fixed effects model as well as various mixed effects
models. The result of the below is one fixed effects model for biomass
from the unique combination of predictor variables. This is repeated for
each unique combination and each model is added to a list of all the
models, with the name reflecting the ‘family’ of the model (basically,
which set of predictors was used) as well as the level of fixed or mixed
effects. These names are important for identifying and evaluating model
families in post-hoc analyses.

``` r
###  unseen in this script is that p is one row of the above printed table, basically one set of predictor variables. 
###  The function loops through each of the rows to create different models.
preds <- perms[p,]
##  get the log-transformed version of the variable
log_preds <- paste0("log", preds)
###  remove missing values
data_clean <- dat |> filter(across(c(log_preds,paste0("log", responseVar)), ~ !is.na(.)))
###  create the model formulas, fixed effects first
fixed_formula <- reformulate(log_preds, response = paste0("log", responseVar))
###  run the model
###  first fixed effects
fixed_model <- lm(fixed_formula, data = data_clean)
##   add the model to the list and name it appropriately
mods <- append(mods,list(fixed_model))
names(mods)[length(mods)] <- paste0(model,".Fixed")
###  get summary stats
fixed_summary <- summary(fixed_model)
##  Akaike and Bayesian Information Criteria
AIC_fixed <- AIC(fixed_model)
BIC_fixed <- BIC(fixed_model)
##  Root Mean Squared Error
RMSE_fixed <- sqrt(mean(resid(fixed_model)^2))
## r-squared
R2_fixed <- fixed_summary$r.squared
##  sigma
sig_fixed <-  sigma(fixed_model)
##  variance inflation factor, only valid on multiple regression models
if (k>1){VIF_fixed <- car::vif(fixed_model)}else{VIF_fixed <- NA}
###  predict on the model
predsF <- predsF |> left_join(data.frame(ID=data_clean$ID,pred=predict(fixed_model,newdata=data_clean)),by=join_by(ID))
###  get the coefficient names
nmsf <- c(names(coef(fixed_model)),names(coef(fixed_model))[-1])
###  make the coefficent names more generic so they can be added to a table with other model coefficients
nmsf[2:(2*length(preds)+1)] <- c(paste0("slope.var",1:(length(preds))),paste0("VIF.var",1:(length(preds))))
nmsm <- nmsf
###  create a data frame row to hold the coefficients
###  this will be added to a growing dataframe of model coefficients
fixed_cols <- as_tibble_row(setNames(as.list(c(coef(fixed_model),VIF_fixed)),
                                     paste0(nmsf, "_Fixed")))
```

The result of the above is one row of a dataframe with some standard
information about the model. As an example, here is the result for a
two-variable model of biomass as a function of DBH and height:

    ## # A tibble: 1 × 5
    ##   `(Intercept)_Fixed` slope.var1_Fixed slope.var2_Fixed VIF.var1_Fixed VIF.var2_Fixed
    ##                 <dbl>            <dbl>            <dbl>          <dbl>          <dbl>
    ## 1                1.90            0.618             1.85           5.05           5.05

Next are the more complicated mixed effects models. For every fixed
effects model there are six mixed effects model, three for random
intercepts only (one for each of the grouping variables and another that
has both grouping variables) and three also for random intercepts and
slopes. In the below code, these model structures are created but the
models are not yet run.

``` r
###  mixed effects
###  random intercept, with all grouping variables
mixed_formula_int <- as.formula(paste0(
  deparse(fixed_formula,width.cutoff = 100L),
  " + ",
  paste0("(1|", groupVars, ")", collapse = " + "))
)
###  random intercept and slope, with all grouping variables
mixed_formula_int_slope <- as.formula(paste0(
  deparse(fixed_formula,width.cutoff = 100L),
  " + ",
  paste0("(",log_preds,"|", rep(groupVars,length(preds)), ")", collapse = " + "))
)
###  now add in random effects on individual grouping variables
for (m in 1:length(groupVars)){
  mixed_formula_int <- append(mixed_formula_int,
                              as.formula(paste0(
                                deparse(fixed_formula,width.cutoff = 100L),
                                " + ",
                                "(1|", groupVars[m], ")")))
names(mixed_formula_int)[length(mixed_formula_int)] <- paste0("MixedInt_",groupVars[m])
mixed_formula_int_slope <- append(mixed_formula_int_slope,
                                  as.formula(paste0(
                                    deparse(fixed_formula,width.cutoff = 100L),
                                    " + ",
                                    paste0("(",log_preds,"|", rep(groupVars[m],length(preds)), ")", collapse = " + "))))
names(mixed_formula_int_slope)[length(mixed_formula_int_slope)] <- paste0("MixedIntSlope_",groupVars[m])
}
names(mixed_formula_int)[1] <- paste0("MixedInt_",paste(groupVars,collapse="&"))
names(mixed_formula_int_slope)[1] <- paste0("MixedIntSlope_",paste(groupVars,collapse="&"))
###  create the random effects combination options
mmmods <- c(paste0(paste(groupvars,collapse=""),c("int","intslope")),outer(groupvars,c("int","intslope"),paste0))
```

As an example of the above product, below are all of the model
structures for the first set of models, which use DBH as a predictor.
The firs model, the fixed effects model has already been run at this
point, but the remaining models have not, we have only established their
structure at this point. We will run each one individually in the next
step, as well as collect coefficient and performance information.

    ##                                                                model                         name                                               description
    ## 1                                              logAGB.kg ~ logDBH.cm                      1.Fixed                                       fixed effects model
    ## 2                 logAGB.kg ~ logDBH.cm + (1 | Species) + (1 | Site)      1.MixedInt_Species&Site            both species and location as random intercepts
    ## 3 logAGB.kg ~ logDBH.cm + (logDBH.cm | Species) + (logDBH.cm | Site) 1.MixedIntSlope_Species&Site both species and location as random intercepts and slopes
    ## 4                              logAGB.kg ~ logDBH.cm + (1 | Species)           1.MixedInt_Species                              species as random intercepts
    ## 5                      logAGB.kg ~ logDBH.cm + (logDBH.cm | Species)      1.MixedIntSlope_Species                   species as random intercepts and slopes
    ## 6                                 logAGB.kg ~ logDBH.cm + (1 | Site)              1.MixedInt_Site                             location as random intercepts
    ## 7                         logAGB.kg ~ logDBH.cm + (logDBH.cm | Site)         1.MixedIntSlope_Site                  location as random intercepts and slopes

In the penultimate step, we iterate through the above created mixed
effects model structures and attempt to fit them. We use a ‘tryCatch’
statement to avoid stopping the routine in the case of a model that, for
whatever reason, throws an error. This is more common in mixed effects
models due to their complex structure and requirements. If the fit is
successful, we extract the desired information, including the same
metrics we gathered for the fixed effects version, as well as the
Intraclass Correlation Coefficient, which is important in evaluating
mixed effects models [Snijders & Bosker,
2012](https://www.stats.ox.ac.uk/~snijders/mlbook.htm). Most of the
other metrics are equivalent to the fixed effects metrics, however the
r-squared returned here is the ‘conditional r-squared’. This is the
r-squared value considering both fixed and mixed effects. See
MuMIn::r.squaredGLMM for more information. In the end, we add together
the collected metrics and coefficients for the fixed effects model and
all of the mixed effects model for the current predictor set.

``` r
###  for each of the ME model formula templates, run the combinations
for(mm in 1:length(mixed_formula_int)){
  MMods <- c(mixed_formula_int[mm],mixed_formula_int_slope[mm])
  if (mm==1){ME = paste(groupVars,collapse="&")}else{ME = groupVars[mm-1]}
  ##  Fit mixed model
  ##  create lists to hold coefficients and model metrics
  R2_mixed <- sigs_mixed <- AIC_mixed <-BIC_mixed<-RMSE_mixed<-coefs_mixed <- ICC <- VIF_mixed <- list()
  for (M in 1:length(MMods)){
    ###  fit the model, if possible given the data
    mixed_model <- tryCatch({
      lmer(MMods[M][[1]], data = data_clean, REML = FALSE)
      }, error = function(e) NULL)
    ###  if the model ran successfully, get the coefficients and metrics
    if (!is.null(mixed_model)) {
      ###  first the predictions
      predsMM <- predsMM |> left_join(data.frame(ID=data_clean$ID,pred=predict(mixed_model,newdata=data_clean)),by=join_by(ID))
      ##  add the model to the list and name it appropriately
      mods <- append(mods,list(mixed_model))
      names(mods)[length(mods)] <- paste(model,names(MMods[M]),sep=".")
      # r-squared, in this case we are using the 'conditional r-squared', which includes both fixed and random effects 
      R2_mixed <- append(R2_mixed,MuMIn::r.squaredGLMM(mixed_model)[, "R2c"])
      # sigma
      sigs_mixed <- append(sigs_mixed,sigma(mixed_model))
      # Aikaike and Bayes Information Criteria
      AIC_mixed <- append(AIC_mixed,AIC(mixed_model))
      BIC_mixed <- append(BIC_mixed,BIC(mixed_model))
      ## Root mean squared error
      RMSE_mixed <- append(RMSE_mixed ,sqrt(mean(resid(mixed_model)^2)))
      #  variance inflation factor
      if (k>1){VIF_mixed <- car::vif(mixed_model)}else{VIF_mixed <-NA}
      ##  add the coefficients together
      coefs_mixed <- append(coefs_mixed,c(fixef(mixed_model,add.dropped=TRUE),VIF_mixed))
      #  get the intraclass correlation coefficient
      ICC <- append(ICC,icc(mixed_model)[1])
      } else {
      ####  if the model could not run, set its metrics and coefficients to NA
      predsMM <- cbind(predsMM,rep(NA,nrow(predsMM)))
      R2_mixed <- append(R2_mixed,NA)
      sigs_mixed <- append(sigs_mixed,NA)
      AIC_mixed <- append(AIC_mixed,NA)
      BIC_mixed <- append(BIC_mixed,NA)
      RMSE_mixed <- append(RMSE_mixed,NA)
      ICC <- append(ICC,NA)
      VIF_mixed <- append(VIF_mixed,NA)
      coefs_mixed[[M]] <- setNames(rep(NA, length(log_preds) + 1),
                                   c("(Intercept)", log_preds))
      }
  }
  ### with both the fixed and mixed effects models fit for the same set of predictors, we add their collected information to a dataset where
  ##  they can be more easily compared. 
  result_row <- tibble(
    VarGroup = varGroup,
    Model = paste(preds, collapse = ", "),
    NumPredictors = length(preds),
    Rsq_Fixed = R2_fixed,
    Sig_Fixed = sig_fixed,
    AIC_Fixed = AIC_fixed,
    BIC_Fixed = BIC_fixed,
    RMSE_Fixed = RMSE_fixed,
    ) |> mutate(MixedEffects=ME) |>
    bind_cols(as_tibble_row(setNames(as.list(c(R2_mixed,sigs_mixed,AIC_mixed,BIC_mixed,RMSE_mixed,ICC)),
                                     paste0(rep(c("Rsq_Mixed","Sig_Mixed","AIC_Mixed","BIC_Mixed","RMSE_Mixed","ICC_Mixed"),each=2),
                                            c("Int","IntSlope")))))
  coef_row <- tibble(
    VarGroup = varGroup,
    MixedEffects=ME,
    Model = paste(preds, collapse = ", "),
    NumPredictors = length(preds)) |>
    bind_cols(fixed_cols)|>
    bind_cols(as_tibble_row(setNames(as.list(coefs_mixed),
                                     paste0(nmsm, "_Mixed",rep(c("Int","IntSlope"),each=length(nmsm))))))
  results[[length(results) + 1]] <- result_row
  coefs[[length(coefs) + 1]] <- coef_row
}
```

The result of the above is a dataframe containing relevant metrics and
model coefficents for all of the models within a given ‘family’, that is
all of the models that use the same set of predictors but with different
combinations of mixed effects. This is done for all of the families, and
the results are combined together into a master data frame containing
all of the information from all of the models.

    ## # A tibble: 3 × 21
    ##   VarGroup Model            NumPredictors Rsq_Fixed Sig_Fixed AIC_Fixed BIC_Fixed RMSE_Fixed MixedEffects Rsq_MixedInt Rsq_MixedIntSlope Sig_MixedInt Sig_MixedIntSlope AIC_MixedInt AIC_MixedIntSlope BIC_MixedInt BIC_MixedIntSlope RMSE_MixedInt RMSE_MixedIntSlope ICC_MixedInt ICC_MixedIntSlope
    ##      <dbl> <chr>                    <dbl>     <dbl>     <dbl>     <dbl>     <dbl>      <dbl> <chr>               <dbl>             <dbl>        <dbl>             <dbl>        <dbl>             <dbl>        <dbl>             <dbl>         <dbl>              <dbl>        <dbl> <lgl>            
    ## 1        6 Height.m, DBH.cm             2     0.955     0.458      523.      539.      0.457 Species&Site        0.966             0.979        0.401             0.322         442.              319.         466.              383.         0.396              0.316        0.229 NA               
    ## 2        6 Height.m, DBH.cm             2     0.955     0.458      523.      539.      0.457 Species             0.966             0.979        0.401             0.322         442.              319.         466.              383.         0.396              0.316        0.229 NA               
    ## 3        6 Height.m, DBH.cm             2     0.955     0.458      523.      539.      0.457 Site                0.966             0.979        0.401             0.322         442.              319.         466.              383.         0.396              0.316        0.229 NA

    ## # A tibble: 3 × 19
    ##   VarGroup MixedEffects Model            NumPredictors `(Intercept)_Fixed` slope.var1_Fixed slope.var2_Fixed VIF.var1_Fixed VIF.var2_Fixed `(Intercept)_MixedInt` slope.var1_MixedInt slope.var2_MixedInt VIF.var1_MixedInt VIF.var2_MixedInt `(Intercept)_MixedIntSlope` slope.var1_MixedIntSlope slope.var2_MixedIntSlope VIF.var1_MixedIntSlope VIF.var2_MixedIntSlope
    ##      <dbl> <chr>        <chr>                    <dbl>               <dbl>            <dbl>            <dbl>          <dbl>          <dbl>                  <dbl>               <dbl>               <dbl>             <dbl>             <dbl>                       <dbl>                    <dbl>                    <dbl>                  <dbl>                  <dbl>
    ## 1        6 Species&Site Height.m, DBH.cm             2                1.90            0.618             1.85           5.05           5.05                   1.87               0.656                1.83              6.04              6.04                        1.92                    0.504                     1.95                   1.17                   1.17
    ## 2        6 Species&Site Height.m, DBH.cm             2                1.90            0.618             1.85           5.05           5.05                   1.87               0.656                1.83              6.04              6.04                        1.92                    0.504                     1.95                   1.17                   1.17
    ## 3        6 Species&Site Height.m, DBH.cm             2                1.90            0.618             1.85           5.05           5.05                   1.87               0.656                1.83              6.04              6.04                        1.92                    0.504                     1.95                   1.17                   1.17
