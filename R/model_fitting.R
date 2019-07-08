#    Copyright (c) 2019 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
#  
#    This file is part of the rCPDMS program.
#
#    rCPDMS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



## Model fitting functions

utils::globalVariables(c("."))

# Pulldown ----------------------------------------------------------------

fitAll_pulldown = function(dat_input, baseline) {
  # dat_input is either a list of dataframes or a single dataframe
  if (!is.data.frame(dat_input)) {
    start_time <- Sys.time()
    cat("================= Fitting models for all the proteins ================= \n")
    results_list = lapply(dat_input, fitAll_pulldown, baseline)
    end_time <- Sys.time()
    print(end_time - start_time)
    return(results_list)
  } else {
    this_project = unique(as.character(dat_input$Project_Experiment_SampleGroup))
    cat("+++ Project_Experiment_SampleGroup: ",
        this_project,
        " +++ \n")
    dat_split = split(dat_input, dat_input$proteinID)
    dat_split = dat_split[lapply(dat_split, nrow) > 0]
    results_all_list = lapply(dat_split, fit_one_protein_pulldown, baseline)
    results_all = bind_rows(results_all_list)
    # significant B: if less than 4 doses
    results_all = sigB_pulldown(dat_input, results_all, baseline)
    # Limma:
    limma_out <- prolimma(dat_input, 'Dose_F', baseline)
    results_all = merge(results_all, limma_out, by="proteinID", all=TRUE)
    # add pVal_min, pVal_limma_min, pVal_sigB_min
    pval_cols = names(results_all)[grepl("pVal_", names(results_all))]
    results_all$pVal_min = apply(results_all[pval_cols], 1, min, na.rm = TRUE)
    pval_cols = names(results_all)[grepl("pVal_limma_", names(results_all))]
    results_all$pVal_limma_min = apply(results_all[pval_cols], 1, min, na.rm = TRUE)
    pval_cols = names(results_all)[grepl("pVal_sigB_", names(results_all))]
    results_all$pVal_sigB_min = apply(results_all[pval_cols], 1, min, na.rm = TRUE)
    return(results_all)
  }
}

fit_one_protein_pulldown = function(dat_one_protein, baseline) {
  op <- options(show.error.messages = FALSE, warn = -1)
  on.exit(options(op))
  
  # dat_one_protein if no missing, number of rows equals to the number of SampleName
  this_protein_id = as.character(unique(dat_one_protein$proteinID))
  
  # results for one protein will be a one-row dataframe
  if ("raw_log_Result_Value" %in% colnames(dat_one_protein)) {
    results = as_tibble(dat_one_protein) %>%
      group_by(Project_Experiment_SampleGroup, proteinID) %>%
      summarise(count = sum(!is.na(log_Result_Value)), count_raw = sum(!is.na(raw_log_Result_Value)))
  } else {
    results = as_tibble(dat_one_protein) %>%
      group_by(Project_Experiment_SampleGroup, proteinID) %>%
      summarise(count = sum(!is.na(log_Result_Value)))
  }
  dat_one_protein = as_tibble(dat_one_protein) %>% filter(!is.na(log_Result_Value))
  if(nrow(dat_one_protein)==0){
    results$note = "Data not available."
    return(results)
  }
  
  remaining_dose = sort(unique(dat_one_protein$Dose))
  if (!(baseline %in% remaining_dose)) {
    cat("Skip fitting ",
        this_protein_id,
        " due to missing baseline dose.\n")
    results$note = "Missing baseline dose."
    return(results)
  } else if (!(length(remaining_dose) > 1)) {
    cat("Skip fitting ",
        this_protein_id,
        " since it does not have more than one dose.\n")
    results$note = "Missing non-baseline dose."
    return(results)
  }
  
  # Model 1: linear regression, dose as numeric log(Dose)
  #           fit when at least 3 doses remaining
  if(length(remaining_dose) >= 3){
    model_1 = try(lm(log_Result_Value ~ logDose, dat_one_protein))
    if (!inherits(model_1, 'try-error')) {
      c_1 = summary(model_1)$coefficients
      results[, c("slope", "pVal_lm")] = c_1[rownames(c_1) == "logDose", c("Estimate", "Pr(>|t|)")]
    }
  }
  
  # Model 2: oneway ANOVA, dose as factor
  model_2 = try(lm(log_Result_Value ~ Dose_F, dat_one_protein))
  if (!inherits(model_2, 'try-error')) {
    c_2 = summary(model_2)$coefficients
    for (d in unique(dat_one_protein$Dose_F[dat_one_protein$Dose > 0])) {
      results[, paste0("est_logRatioDose_", d)] = c_2[paste0("Dose_F", d), "Estimate"]
      results[, paste0("pVal_logRatioDose_", d)] = c_2[paste0("Dose_F", d), "Pr(>|t|)"]
    }
    results$pVal_ANOVA =  anova(model_2)$"Pr(>F)"[1]
  }
  
  
  # Model 3: logistic curve
  #           fit when more than 4 data points and at least 4 doses present
  if (results$count > 4 & length(remaining_dose) >= 4) {
    model_3 = try(drc::drm(log_Result_Value ~ logDose, data = dat_one_protein,
      # four-parameter log-logistic models
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
      control = drmc(errorm = FALSE,useD = TRUE,warnVal = 0,rmNA = TRUE,method = "Nelder-Mead")
    ))
    
    if (!inherits(model_3, 'try-error') &&
        class(model_3) == "drc" && is.null(model_3$convergence)) {
      sm = summary(model_3)$coefficients
      results$pVal_logistic = noEffect(model_3)[['p-value']]
      results[, c("logistic_IC50_estimate", "logistic_IC50_Std")] = sm["ED50:(Intercept)", c("Estimate", "Std. Error")]
      results$logistic_slope = sm["Slope:(Intercept)", "Estimate"]
      results$logistic_lowerLimit = sm["Lower Limit:(Intercept)", "Estimate"]
      results$logistic_upperLimit = sm["Upper Limit:(Intercept)", "Estimate"]
      # quality control of estimates
      if (!is.na(results[,"logistic_IC50_estimate"])){
        if(results[,"logistic_IC50_estimate"] > max(dat_one_protein$logDose)){
          results[, "logistic_IC50_estimate"] = max(dat_one_protein$logDose)
          results[, "note"] = "logistic_IC50_estimate > max logDose"
        }
        if(results[,"logistic_IC50_estimate"] < min(dat_one_protein$logDose)){
          results[, "logistic_IC50_estimate"] = min(dat_one_protein$logDose)
          results[, "note"] = "logistic_IC50_estimate < min logDose"
        }
      }
    }
  }else{
    # add empty columns if logistic model is not fitted.
    results$pVal_logistic = NA
    results$logistic_IC50_estimate = NA
    results$logistic_IC50_Std = NA
    results$logistic_slope = NA
    results$logistic_lowerLimit = NA
    results$logistic_upperLimit = NA
  }
  
  # add maxLogRatio
  #results$maxLogRatio = max(dat_one_protein$log_Result_Value)-min(dat_one_protein$log_Result_Value)
  # revise maxLogRatio to show the largest est effect
  x = as.numeric(results[, which(startsWith(names(results), "est_logRatioDose_"))])
  if (length(x) > 0) {
    #results$maxLogRatio = x[max.col(abs(x))]
    results$maxLogRatio = x[which.max(abs(x))]
  } else {
    results$maxLogRatio = NA
  }
  
  return(results)
}

sigB_pulldown = function(dat_input, result, baseline, DoseName = "Dose"){
  dat_input$Dose = dat_input[[DoseName]]
  # run fitting only when no more than 3 doses
  if(length(unique(dat_input$Dose)) >3 ){return(result)}
  nonBaselineDoses = setdiff(unique(dat_input[[DoseName]]),baseline)
  # if design is single repeat, run significant B method, otherwise initialize NA columns
  if(length(unique(dat_input$Repeat[!is.na(dat_input$Repeat)]))==1){
    for(d in nonBaselineDoses){
      # first binning the proteins according to avg intensity between DMSO and treatment
      temp = dat_input %>% filter(Dose %in% c(d, baseline)) %>% 
        select(proteinID, log_Result_Value) %>%
        group_by(proteinID) %>% summarise(avg = mean(log_Result_Value, na.rm = FALSE)) %>% 
        ungroup() %>% filter(!is.na(avg)) %>% arrange(avg)
      n_bins = floor(nrow(temp)/300)
      if(n_bins>=1){
        n_per_bins = ceiling(nrow(temp)/n_bins)
        temp$group = rep(1:n_bins, each = n_per_bins)[1:nrow(temp)]
      }else{
        temp$group = 1
      }
      temp$avg = NULL
      temp = temp %>% inner_join(result[,c("proteinID", paste0("est_logRatioDose_",d))], by = "proteinID")
      temp = significance_B(temp)
      result = result %>% left_join(temp[,c("proteinID", paste0("pVal_sigB_logRatioDose_",d))], by = "proteinID")
    }
  }else{
    logRatio_cols_inds = grep("est_logRatioDose_", names(result))
    newCols = gsub("est_","pVal_sigB_",names(result)[logRatio_cols_inds])
    result[newCols] = NA
  }
  return(result)
}

prolimma <- function(dat, designVar, baseline) {
  dat$TrtRep <- paste(dat[[designVar]], dat$Repeat, sep=":")
  dat_wide <- spread(dat[, c('proteinID', 'TrtRep', 'log_Result_Value')], 
                     TrtRep, log_Result_Value)
  Trt <- factor(sapply(strsplit(names(dat_wide)[-1], ':'), '[', 1))
  origLvls <- levels(Trt) 
  levels(Trt) <- make.names(levels(Trt))
  Trts <- levels(Trt)
  baseline <- make.names(baseline)
  baseLevel <- which(Trts == baseline)
  trtLevel <- which(Trts != baseline)
  design <- model.matrix(~ 0 + Trt)
  colnames(design) <- Trts
  cntrst <- as.list(paste(Trts[trtLevel], '-', Trts[baseLevel]))
  cmat <- do.call(makeContrasts, c(cntrst, list(levels=design)))
  fit <- lmFit(dat_wide[, -1], design)
  fit2 <- contrasts.fit(fit, cmat)
  fit2 <- eBayes(fit2)
  p_limma <- fit2$p.value
  lvl <- sapply(strsplit(colnames(p_limma), ' '), '[', 1)
  colnames(p_limma) <- paste('pVal_limma_', origLvls[match(lvl, levels(Trt))], sep='')
  cbind(data.frame(proteinID=dat_wide$proteinID), p_limma)
}

# Pro10K ------------------------------------------------------------------

fitAll_pro10K = function(dat_input, baseline) {
  # dat_input is either a list of dataframes or a single dataframe
  if (!is.data.frame(dat_input)) {
    start_time <- Sys.time()
    cat("================= Fitting models for all the proteins ================= \n")
    results_list = lapply(dat_input, fitAll_pro10K, baseline)
    end_time <- Sys.time()
    print(end_time - start_time)
    return(results_list)
  } else {
    this_project = unique(as.character(dat_input$Project_Experiment))
    cat("+++ Project_Experiment: ", this_project, " +++ \n")
    
    dat_split = split(dat_input, dat_input$proteinID)
    dat_split = dat_split[lapply(dat_split, nrow) > 0]
    results_all_list = lapply(dat_split, fit_one_protein_pro10K, baseline)
    results_all = bind_rows(results_all_list)
    
    # Limma:
    limma_p <- prolimma(dat_input, 'SampleGroup', baseline)
    limma_p <- gather(limma_p, SampleGroup, pVal_limma, -proteinID)
    limma_p$SampleGroup <- sub('pVal_limma_', '', limma_p$SampleGroup)
    results_all <- merge(results_all, limma_p, by=c('proteinID', 'SampleGroup'), all=TRUE)
    # add pVal_min for each row (protein*SampleGroup)
    pval_cols = names(results_all)[grepl("pVal_", names(results_all))]
    results_all$pVal_min = apply(results_all[pval_cols], 1, min, na.rm = TRUE)
    results_all$pVal_min[results_all$SampleGroup=="DMSO"] = NA 
    return(results_all)
  }
}

fit_one_protein_pro10K = function(dat_one_protein, baseline = "DMSO") {
  op <- options(show.error.messages = FALSE, warn=-1)
  on.exit(options(op))
  # dat_one_protein if no missing, number of rows equals to the number of wells
  this_protein_id = as.character(unique(dat_one_protein$proteinID))
  dat_one_protein = as_tibble(dat_one_protein) %>% filter(!is.na(log_Result_Value))
  
  if ("raw_log_Result_Value" %in% colnames(dat_one_protein)) {
    results = as_tibble(dat_one_protein) %>%
      group_by(Project_Experiment, proteinID, SampleGroup) %>%
      summarise(count = n(), count_raw = sum(!is.na(raw_log_Result_Value))) %>%
      mutate(SampleGroup = as.character(SampleGroup))
  } else{
    results = as_tibble(dat_one_protein) %>%
      group_by(Project_Experiment, proteinID, SampleGroup) %>%
      summarise(count = n()) %>%
      mutate(SampleGroup = as.character(SampleGroup))
  }
  
  if (sum(dat_one_protein$SampleGroup == baseline) == 0) {
    cat("No baseline found in protein: ", this_protein_id, "\n")
    #results = data.frame(proteinID = this_protein_id,
    #                     SampleGroup = "DMSO",
    #                     count = 0, log2_Ratio = NA, pVal = NA, stringsAsFactors = FALSE)
    #results = as_tibble(results)
    results = results %>% mutate(log2_Ratio = as.numeric(NA), pVal = as.numeric(NA))
    
  } else if (sum(dat_one_protein$SampleGroup != baseline) == 0) {
    cat("No compounds found in protein: ", this_protein_id, "\n")
    results = results %>% mutate(log2_Ratio = as.numeric(NA), pVal = as.numeric(NA))
  } else {
    model = lm(log_Result_Value ~ SampleGroup, dat_one_protein)
    coefMat = summary(model)$coefficients
    #results = data.frame(proteinID = rep(this_protein_id, nrow(coefMat)))
    df_coefMat = data.frame(
      SampleGroup = sapply(rownames(coefMat), function(x) gsub("SampleGroup", "", x)),
      stringsAsFactors = FALSE
    )
    #results$drug = sapply(rownames(coefMat), function(x) gsub("SampleGroup", "", x))
    #results$drug[which(results$drug=="(Intercept)")] = "DMSO"
    df_coefMat$SampleGroup[which(df_coefMat$SampleGroup == "(Intercept)")] = baseline
    df_coefMat$log2_Ratio = coefMat[, "Estimate"] #* log2(10)
    df_coefMat$pVal = coefMat[, "Pr(>|t|)"]
    results = results %>% left_join(df_coefMat, by = "SampleGroup")
  }
  return(results)
}


# CETSA -------------------------------------------------------------------

fitAll_CETSA = function(dat_input, baseline){
  # dat_input is either a list of dataframes or a single dataframe
  if(!is.data.frame(dat_input)){
    start_time <- Sys.time()
    cat("================= Fitting models for all the proteins ================= \n")
    results_list = lapply(dat_input, fitAll_CETSA, baseline)
    # add IC50 diff and pvalue to compare other groups with DMSO groups
    results_list = compareGroups(results_list, baseline)
    end_time <- Sys.time()
    print(end_time - start_time)
    return(results_list)
  }else{
    this_project = unique(as.character(dat_input$Project_Experiment))
    cat("+++ Project_Experiment: ",this_project," +++ \n")
    
    dat_split = split(dat_input, dat_input$proteinID)
    dat_split = dat_split[lapply(dat_split, nrow) > 0]
    results_all_list = lapply(dat_split, fit_one_protein_CETSA, printWarning = FALSE)
    results_all = bind_rows(results_all_list)
    return(results_all)
  }
}

# fit models for each protein and all samplegroups
fit_one_protein_CETSA = function(dat_one_protein, printWarning = FALSE){
  op <- options(show.error.messages = FALSE, warn=-1)
  on.exit(options(op))
  
  # dat_one_protein if no missing, number of rows equals to #SampleName * #Temperature
  this_protein_id = as.character(unique(dat_one_protein$proteinID))
  dat_one_protein = as_tibble(dat_one_protein) %>% filter(!is.na(log_Result_Value))
  
  if ("log_Result_Value_raw" %in% colnames(dat_one_protein)) {
    results = as_tibble(dat_one_protein) %>% 
      group_by(Project_Experiment, proteinID, SampleGroup) %>% 
      summarise(count = n(), count_raw = sum(!is.na(log_Result_Value_raw))) %>% 
      mutate(SampleGroup = as.character(SampleGroup))
  } else {
    results = as_tibble(dat_one_protein) %>% 
      group_by(Project_Experiment, proteinID, SampleGroup) %>% 
      summarise(count = n()) %>% 
      mutate(SampleGroup = as.character(SampleGroup))
  }
  
  results$IC50_estimate = NA
  results$IC50_Std = NA
  results$pValue = NA
  results$note = NA
  results$maxLogRatio = NA
  results$slope = NA
  results$lowerLimit = NA
  results$upperLimit = NA
  
  for (i in 1:NROW(results)) {
    
    dat_i = dat_one_protein[(dat_one_protein$SampleGroup==results$SampleGroup[i]),c("log_Result_Value","Temperature")]
    nrow_dat_i = nrow(dat_i)
    
    if(nrow_dat_i == 0) {
      results$note[i] = "missing raw data"
      next
    } else if (nrow_dat_i < 6) {
      results$note[i] = "insufficient data"
      next
    } else {
      fit = try(drc::drm(log_Result_Value ~ Temperature, 
                         data = dat_i, 
                         # four-parameter log-logistic models
                         fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                         control = drmc(errorm = FALSE, useD = TRUE, warnVal = 0,
                                        rmNA = TRUE, method = "Nelder-Mead")))
      if (!inherits(fit, 'try-error') && inherits(fit, "drc") && is.null(fit$convergence)) {
        if (exists("fit") && fit$fit$convergence) {
          results$pValue[i] = noEffect(fit)[['p-value']]
          if(noEffect(fit)[["p-value"]] < 0.01) {
            sm = summary(fit)$coefficients
            fittedValues = fitted(fit)
            results[i,c("IC50_estimate","IC50_Std")] = sm[4,1:2]
            
            # to eliminate inaccurate IC50_estimate due to poor data quality
            if (!is.na(results[i,"IC50_Std"]) && results[i,"IC50_Std"] > 10) {
              results[i,"IC50_estimate"] = NA
              results[i,"note"] = "IC50_Std > 10"
            }
            # to eliminate unreasonable IC50_estimate due to poor data quality
            # if (!is.na(results[i,"IC50_estimate"]) && 
            #     (results[i,"IC50_estimate"] > 100 || results[i,"IC50_estimate"] < 0)) {
            #   results[i, "IC50_estimate"] = NA
            #   results[i, "note"] = "IC50_estimate outside (0,100)"
            # }
            if (!is.na(results[i,"IC50_estimate"])){
              if(results[i,"IC50_estimate"] > 100){
                results[i, "IC50_estimate"] = 100
                results[i, "note"] = "IC50_estimate > 100"
              }
              if(results[i,"IC50_estimate"] < 0){
                results[i, "IC50_estimate"] = 0
                results[i, "note"] = "IC50_estimate < 0"
              }
            }
            results[i, "maxLogRatio"] = max(fittedValues) - min(fittedValues)
            results$slope[i] = sm["Slope:(Intercept)", "Estimate"]
            results$lowerLimit[i] = sm["Lower Limit:(Intercept)", "Estimate"]
            results$upperLimit[i] = sm["Upper Limit:(Intercept)", "Estimate"]
          } else {
            results[i, "note"] = "noEffect"
          }
        } else {
          results[i, "note"] = "Fail to converge"
        }
      } else {
        results[i, "note"] = "Error during model fitting"
      }
      if ( "fit" %in% ls()) { rm(fit) }
    }
  }
  return(results)
}

compareGroups = function(input, baseline = "DMSO") {
  # input is either results_list or a results dataframe
  if (is.data.frame(input)) {
    output = input
    output$IC50_diff = NA
    output$pValue_IC50_diff = NA
    for (i in 1:nrow(input)) {
      results_i = input[i,]
      if (results_i$SampleGroup != baseline) {
        results_baseline = input[which(input$proteinID==results_i$proteinID & input$SampleGroup==baseline),]
        if (nrow(results_baseline) == 0) {
          next
        }
        #d_f = results_i$numOfPoints - 4 # since there are 4 parameters
        IC50_diff_i = results_i$IC50_estimate - results_baseline$IC50_estimate
        sss = sqrt(results_i$IC50_Std^2 + results_baseline$IC50_Std^2)
        pValue_IC50_diff_i = 2 * (1 - pnorm(abs(IC50_diff_i) / sss))
        
        output$IC50_diff[i] = IC50_diff_i
        output$pValue_IC50_diff[i] = pValue_IC50_diff_i
      }
    }
  } else {
    output = lapply(input, compareGroups, baseline = baseline)
  }
  return(output)
}


# 2dCETSA -----------------------------------------------------------------

fitAll_2dCETSA = fitAll_CETSA

fitAll_2dCETSA_by_temp = function(dat_input, baseline) {
  # dat_input is either a list of dataframes or a single dataframe
  if (!is.data.frame(dat_input)) {
    start_time <- Sys.time()
    cat("================= Fitting models for all the proteins ================= \n")
    results_list = lapply(dat_input, fitAll_2dCETSA_by_temp, baseline)
    end_time <- Sys.time()
    print(end_time - start_time)
    return(results_list)
  } else{
    this_project = unique(as.character(dat_input$Project_Experiment))
    cat("+++ Project_Experiment: ", this_project, " +++ \n")

    dat_input = add_logDose(dat_input)
    dat_input = dat_input[!is.na(dat_input$log_Result_Value), ]
    dat_split = split(dat_input,
                      list(dat_input$proteinID, dat_input$Temperature))
    dat_split = dat_split[lapply(dat_split, nrow) > 0]
    logdose_range = c(0, max(dat_input$logDose) + 1)
    results_all_list = lapply(dat_split, fit_one_protein_pulldown2, logdose_range)
    results_all = bind_rows(results_all_list)
    # significance B for each temperature
    dat_input$proteinID = as.character(dat_input$proteinID)
    dat_input$SampleGroup = as.character(dat_input$SampleGroup)
    results_all$proteinID = as.character(results_all$proteinID)
    T_list = unique(dat_input$Temperature)
    results_all_list_B = list()
    for(i in 1:length(T_list)){
      dat_t = dat_input[dat_input$Temperature==T_list[i],]
      result_t = results_all[results_all$Temperature==T_list[i],]
      results_all_list_B[[i]] = sigB_pulldown(dat_t, result_t, baseline, DoseName = "SampleGroup")
    }
    results_all = bind_rows(results_all_list_B)
    # add pVal_min, pVal_limma_min, pVal_sigB_min
    pval_cols = names(results_all)[grepl("pVal_", names(results_all))]
    results_all$pVal_min = apply(results_all[pval_cols], 1, min, na.rm = TRUE)
    pval_cols = names(results_all)[grepl("pVal_limma_", names(results_all))]
    results_all$pVal_limma_min = apply(results_all[pval_cols], 1, min, na.rm = TRUE)
    pval_cols = names(results_all)[grepl("pVal_sigB_", names(results_all))]
    results_all$pVal_sigB_min = apply(results_all[pval_cols], 1, min, na.rm = TRUE)
    return(results_all)
  }
}


fit_one_protein_pulldown2 = function(dat_one_protein, logdose_range) {
  op <- options(show.error.messages = FALSE, warn = -1)
  on.exit(options(op))
  
  # dat_one_protein if no missing, number of rows equals to the number of SampleName
  this_protein_id = as.character(unique(dat_one_protein$proteinID))
  dat_one_protein = as_tibble(dat_one_protein) %>% filter(!is.na(log_Result_Value))
  # results for one protein will be a one-row dataframe
  if ("log_Result_Value_raw" %in% colnames(dat_one_protein)) {
    results = as_tibble(dat_one_protein) %>%
      group_by(Project_Experiment, Temperature, proteinID) %>%
      summarise(count = n(), count_raw = sum(!is.na(log_Result_Value_raw)))
  } else {
    results = as_tibble(dat_one_protein) %>%
      group_by(Project_Experiment, Temperature, proteinID) %>%
      summarise(count = n())
  }
  if (nrow(dat_one_protein) == 0) {
    return(results)
  }
  
  remaining_dose = sort(unique(dat_one_protein$Dose))
  if (min(remaining_dose) != 0) {
    cat("Skip fitting ",
        this_protein_id,
        "due to missing baseline Dose. \n")
    results$note = "Missing baseline dose."
    return(results)
  } else if (length(remaining_dose) == 1) {
    cat("Skip fitting ",
        this_protein_id,
        "since it has only one Dose. \n")
    results$note = "Missing non-baseline dose."
    return(results)
  }
  
  # Model 1: linear regression, dose as numeric log(Dose)
  #           fit when at least 3 doses remaining
  if(length(remaining_dose) >= 3){
    model_1 = try(lm(log_Result_Value ~ logDose, dat_one_protein))
    if (!inherits(model_1, 'try-error')) {
      c_1 = summary(model_1)$coefficients
      results[, c("slope", "pVal_lm")] = c_1[rownames(c_1) == "logDose", c("Estimate", "Pr(>|t|)")]
    }
  }
  
  # Model 2: oneway ANOVA, dose as factor
  model_2 = try(lm(log_Result_Value ~ SampleGroup, dat_one_protein))
  if (!inherits(model_2, 'try-error')) {
    c_2 = summary(model_2)$coefficients
    for (d in unique(dat_one_protein$SampleGroup[dat_one_protein$Dose > 0])) {
      results[, paste0("est_logRatioDose_", d)] = c_2[paste0("SampleGroup", d), "Estimate"]
      results[, paste0("pVal_logRatioDose_", d)] = c_2[paste0("SampleGroup", d), "Pr(>|t|)"]
    }
    results$pVal_ANOVA =  anova(model_2)$"Pr(>F)"[1]
  }
  
  
  # Model 3: logistic curve
  #           fit when more than 4 data points and at least 4 doses present
  if (results$count > 4 & length(remaining_dose) >= 4) {
    model_3 = try(drc::drm(log_Result_Value ~ logDose,data = dat_one_protein,
      # four-parameter log-logistic models
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
      control = drmc(errorm = FALSE,useD = TRUE,warnVal = 0,rmNA = TRUE,method = "Nelder-Mead")
    ))
    if (!inherits(model_3, 'try-error') && inherits(model_3, "drc") && is.null(model_3$convergence)) {
      sm = summary(model_3)$coefficients
      # only keep meaningful results: logistic_IC50_estimate is within a Dose range that make sense
      if ((sm["ED50:(Intercept)", "Estimate"] >= logdose_range[1]) &&
          (sm["ED50:(Intercept)", "Estimate"] <= logdose_range[2])) {
        results$pVal_logistic = noEffect(model_3)[['p-value']]
        results[, c("logistic_IC50_estimate", "logistic_IC50_Std")] = sm["ED50:(Intercept)", c("Estimate", "Std. Error")]
        results$logistic_slope = sm["Slope:(Intercept)", "Estimate"]
        results$logistic_lowerLimit = sm["Lower Limit:(Intercept)", "Estimate"]
        results$logistic_upperLimit = sm["Upper Limit:(Intercept)", "Estimate"]
      }
    }
    # quality control of estimates
    if ("logistic_IC50_estimate" %in% names(results)){
      if(results[,"logistic_IC50_estimate"] > max(dat_one_protein$logDose)){
        results[, "logistic_IC50_estimate"] = max(dat_one_protein$logDose)
        results[, "note"] = "logistic_IC50_estimate > max logDose"
      }
      if(results[,"logistic_IC50_estimate"] < min(dat_one_protein$logDose)){
        results[, "logistic_IC50_estimate"] = min(dat_one_protein$logDose)
        results[, "note"] = "logistic_IC50_estimate < min logDose"
      }
    }
  }else{
    results$pVal_logistic = NA
    results$logistic_IC50_estimate = NA
    results$logistic_IC50_Std = NA
    results$logistic_slope = NA
    results$logistic_lowerLimit = NA
    results$logistic_upperLimit = NA
  }
  
  # add maxLogRatio
  #results$maxLogRatio = max(dat_one_protein$log_Result_Value)-min(dat_one_protein$log_Result_Value)
  x = as.numeric(results[, which(startsWith(names(results), "est_logRatioDose_"))])
  results$maxLogRatio = if (length(x) > 0) x[which.max(abs(x))] else NA
  return(results)
}
