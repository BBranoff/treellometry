#' @importFrom performance icc
#' @importFrom lme4 lmer fixef isSingular VarCorr
#' @importFrom tibble as_tibble_row tibble
#' @importFrom dplyr mutate filter across left_join select n bind_cols join_by bind_rows cur_group_id last_col relocate group_by ungroup
explore_allom_models <- function(dat, responseVar, predictorVars, groupVars,scle=FALSE,varorder=FALSE) {
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
  # Loop through predictor combinations
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
    ###  for each of the predictor sets, run the models
    for (predss in combos) {
      ### if desired, establish all potential combinations of the current predictor set (change the order of predictors).
      if (varorder) perms <- gtools::permutations(length(predss),length(predss),predss)
      ## otherwise, just take the single combination of the variables, regardless of order
      else perms=matrix(predss,ncol=length(predss))
      if (nrow(perms)==0){perms <- matrix(1,nrow = 1,ncol=1)}
      ###  establish the current variable group ID
      varGroup = varGroup + 1
      ##  for each of the predictor combinations, run the models
      for (p in 1:nrow(perms)){
        preds <- perms[p,]
        log_preds <- paste0("log", preds)#,"_scaled")
        ###  remove missing values
        data_clean <- dat |> filter(across(c(log_preds,paste0("log", responseVar)), ~ !is.na(.)))
        ###  create the model formulas, fixed effects first
        fixed_formula <- reformulate(log_preds, response = paste0("log", responseVar))
        ###  run the model
        ###  first fixed effects
        fixed_model <- lm(fixed_formula, data = data_clean)
        mods <- append(mods,list(fixed_model))
        model <- ModelME <- model+1
        names(mods)[length(mods)] <- paste0(model,".Fixed")
        ###  get summary stats
        fixed_summary <- summary(fixed_model)
        cv <- cv::cv(fixed_model, data_clean, k=10, details=TRUE)
        #cv <- performance::performance_cv(eval(substitute(lm(fixed_model,data=data_clean))),method="k_fold",k=10)
        AIC_fixed <- AIC(fixed_model)
        BIC_fixed <- BIC(fixed_model)
        RMSE_fixed <- sqrt(mean(resid(fixed_model)^2))
        RMSE_fixed_std <-  RMSE_fixed/mean(data_clean[,paste0("log", responseVar)],na.rm=TRUE)
        RMSE_CVmean_fixed <-  mean(sqrt(cv$details$criterion),na.rm=TRUE)
        RMSE_CVsd_fixed <-  sd(sqrt(cv$details$criterion),na.rm=TRUE)
        R2_fixed <- fixed_summary$r.squared
        sig_fixed <-  sigma(fixed_model)
        if (k>1){VIF_fixed <- car::vif(fixed_model)}else{VIF_fixed <- NA}
        ###  predict on the model
        predsF <- predsF |> left_join(data.frame(ID=data_clean$ID,pred=predict(fixed_model,newdata=data_clean)),by=join_by(ID))
        ###  get the coefficient names
        nmsf <- c(names(coef(fixed_model)),names(coef(fixed_model))[-1])
        ###  make the coefficent names more generic
        nmsf[2:(2*length(preds)+1)] <- c(paste0("slope.var",1:(length(preds))),paste0("VIF.var",1:(length(preds))))
        nmsm <- nmsf
        ###  create a data frame row to hold the coefficients
        ###  this will be added to a growing dataframe of model coefficients
        fixed_cols <- as_tibble_row(setNames(as.list(c(coef(fixed_model),VIF_fixed)),
                                             paste0(nmsf, "_Fixed")))
        ###  then mixed effects
        ###  random intercept, with all grouping variables
        mixed_formula_int <- as.formula(paste0(
          deparse(fixed_formula,width.cutoff = 100L),
          " + ",
          paste0("(1|", groupVars, ")", collapse = " + ")
        ))
        ###  random intercept and slope, with all grouping variables
        mixed_formula_int_slope <- as.formula(paste0(
          deparse(fixed_formula,width.cutoff = 100L),
          " + ",
          paste0("(",log_preds,"|", rep(groupVars,length(preds)), ")", collapse = " + ")
        ))
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
        ###  now the mixed effects models
        ###  create the random effects combination options
        mmmods <- c(paste0(paste(groupvars,collapse=""),c("int","intslope")),outer(groupvars,c("int","intslope"),paste0))
        ###  for each of the ME model formula templates, run the combinations
        for(mm in 1:length(mixed_formula_int)){
          MMods <- c(mixed_formula_int[mm],mixed_formula_int_slope[mm])
          if (mm==1){ME = paste(groupVars,collapse="&")}else{ME = groupVars[mm-1]}
          ##  Fit mixed model
          ##  create lists to hold coefficients and model metrics
          R2_mixed <- sigs_mixed <- AIC_mixed <-BIC_mixed<-RMSE_mixed<-RMSE_mixed_std<-RMSE_CVmean_mixed<-RMSE_CVsd_mixed<-coefs_mixed <- ICC<-ICC2 <- VIF_mixed <- sings_mixed <- list()
          for (M in 1:length(MMods)){
            ###  fit the model, if possible given the data
            #ModelME <- ModelME+0.1
            mixed_model <- tryCatch({
              lmer(MMods[M][[1]], data = data_clean, REML = FALSE)
            }, error = function(e) NULL)
            ###  if the model ran successfully, get the coefficients and metrics
            if (!is.null(mixed_model)) {
              ##  check for singularity
              sings_mixed <- append(sings_mixed,isSingular(mixed_model))
              ###  first the predictions
              predsMM <- predsMM |> left_join(data.frame(ID=data_clean$ID,pred=predict(mixed_model,newdata=data_clean)),by=join_by(ID))
              ####  caret offers more control but is not optimal for me models
              cv <- cv::cv(mixed_model, data_clean, k=10, details=TRUE)
              mods <- append(mods,list(mixed_model))
              names(mods)[length(mods)] <- paste(model,names(MMods[M]),sep=".")
              R2_mixed <- append(R2_mixed,MuMIn::r.squaredGLMM(mixed_model)[, "R2c"])
              sigs_mixed <- append(sigs_mixed,sigma(mixed_model))
              AIC_mixed <- append(AIC_mixed,AIC(mixed_model))
              BIC_mixed <- append(BIC_mixed,BIC(mixed_model))
              RMSE_mixed <- append(RMSE_mixed ,sqrt(mean(resid(mixed_model)^2)))
              RMSE_mixed_std = append(RMSE_mixed_std,sqrt(mean(resid(mixed_model)^2))/mean(data_clean[, paste0("log", responseVar)],na.rm=TRUE))
              RMSE_CVmean_mixed <- append(RMSE_CVmean_mixed , mean(sqrt(cv$details$criterion),na.rm=TRUE))
              RMSE_CVsd_mixed <- append(RMSE_CVsd_mixed , sd(sqrt(cv$details$criterion),na.rm=TRUE))
              if (k>1){VIF_mixed <- car::vif(mixed_model)}else{VIF_mixed <-NA}
              coefs_mixed <- append(coefs_mixed,c(fixef(mixed_model,add.dropped=TRUE),VIF_mixed))
              #if (is.na(icc(mixed_model)[1])){browser()}
              ICC <- append(ICC,icc(mixed_model)[1])
              vars <- as.data.frame(VarCorr(mixed_model));group_var <- sum(vars[vars$grp != "Residual", "vcov"]);resid_var <- vars[vars$grp == "Residual", "vcov"]
              ICC2 <- append(ICC2,(group_var / (group_var + resid_var)))
            } else {
              ####  if the model could not run, set its metrics and coefficients to NA
              sings_mixed <- append(sings_mixed,NA)
              predsMM <- cbind(predsMM,rep(NA,nrow(predsMM)))
              R2_mixed <- append(R2_mixed,NA)
              sigs_mixed <- append(sigs_mixed,NA)
              AIC_mixed <- append(AIC_mixed,NA)
              BIC_mixed <- append(BIC_mixed,NA)
              RMSE_mixed <- append(RMSE_mixed,NA)
              RMSE_mixed_std <- append(RMSE_mixed_std,NA)
              RMSE_CVmean_mixed <- append(RMSE_CVmean_mixed,NA)
              RMSE_CVsd_mixed <- append(RMSE_CVsd_mixed,NA)
              ICC <- append(ICC,NA)
              ICC2 <- append(ICC2,NA)
              VIF_mixed <- append(VIF_mixed,NA)
              coefs_mixed[[M]] <- setNames(rep(NA, length(log_preds) + 1),
                                           c("(Intercept)", log_preds))
            }
          }
            result_row <- tibble(
            VarGroup = varGroup,
            Model = paste(preds, collapse = ", "),
            NumPredictors = length(preds),
            Rsq_Fixed = R2_fixed,
            Sig_Fixed = sig_fixed,
            AIC_Fixed = AIC_fixed,
            BIC_Fixed = BIC_fixed,
            RMSE_Fixed = RMSE_fixed,
            RMSE.Std_Fixed = RMSE_fixed_std,
            RMSE.CVmean_Fixed = RMSE_CVmean_fixed,
            RMSE.CVsd_Fixed = RMSE_CVsd_fixed,
          ) |> mutate(MixedEffects=ME) |>
            bind_cols(as_tibble_row(setNames(as.list(c(sings_mixed,R2_mixed,sigs_mixed,AIC_mixed,BIC_mixed,RMSE_mixed,RMSE_mixed_std,RMSE_CVmean_mixed,RMSE_CVsd_mixed,ICC,ICC2)),
                                             paste0(rep(c("Singular_Mixed","Rsq_Mixed","Sig_Mixed","AIC_Mixed","BIC_Mixed","RMSE_Mixed","RMSE.Std_Mixed","RMSE.CVmean_Mixed","RMSE.CVsd_Mixed","ICC_Mixed","ICC2_Mixed"),each=2),
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
      }
    }
  }
  #modnamesF <- paste0("mod",(ncol(data_clean)+1):ncol(predsF)-ncol(data_clean),"_")
  #modnamesMM <- paste0("mod",rep((ncol(data_clean)+1):ncol(predsF)-ncol(data_clean),each=6),"_",mmmods)
  names(predsF)[(ncol(data_clean)+1):ncol(predsF)] <- paste0("preds_",grep("Fixed",names(mods),value=TRUE))
  names(predsMM)[(ncol(data_clean)+1):ncol(predsMM)] <- paste0("preds_",grep("Fixed",names(mods),value=TRUE,invert=TRUE))
 # names(mods) <- lapply(seq(1,length(modnamesF)),function(i,nmsF,nmsMM){
   # mmseq <- seq(1,length(nmsMM),by=6)
   # namesnew <- c(nmsF[i],nmsMM[mmseq[i]:(mmseq[i]+5)])
   # namesnew
  #},
 # nmsF=modnamesF,nmsMM=modnamesMM)|>unlist()
  results_df <- bind_rows(results) |>
    mutate(Model=factor(Model,levels=unique(Model))) |>
    group_by(Model)|>
    mutate(ModelN = cur_group_id())|>
    ungroup()|>
    relocate(any_of(contains("Rsq")),.after=NumPredictors)|>
    relocate(any_of(contains("AIC")),.after=Rsq_MixedIntSlope)|>
    relocate(any_of(contains("BIC")),.after=AIC_MixedIntSlope)|>
    relocate(any_of(contains("RMSE")),.after=BIC_MixedIntSlope) |>
    relocate(any_of(contains("MixedEffects")),.after=NumPredictors) |>
    relocate(any_of(contains("ModelN")),.after=Model)

  coef_df <- bind_rows(coefs) |>
    mutate(Model=factor(Model,levels=unique(Model))) |>
    group_by(Model)|>
    mutate(ModelN = cur_group_id()) |>
    ungroup() |>
    relocate(any_of(contains("Intercept")),.after=NumPredictors) |>
    relocate(any_of(contains("VIF")),.after=last_col())|>
    relocate(any_of(contains("ModelN")),.after=Model)
  return(list(results=results_df,coefs=coef_df,predsF=predsF,predsMM=predsMM,mods=mods))
}


