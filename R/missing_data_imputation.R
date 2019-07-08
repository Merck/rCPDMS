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



## Missing data imputation

utils::globalVariables(c(".","log_Result_Value","proteinID","protein_median","SampleName",
                         "raw_log_Result_Value","outlier","NumProteins","Group","..density.."))

# Missing data imputation
#   call function impute_oneProject, which maybe different for different workflow
impute_list = function(workflow, dat_list, saveFolder) {
  start_time <- Sys.time()
  # re-define impute_oneProject function according to workflow
  #impute_oneProject = eval(paste0("impute_oneProject_", workflow))
  impute_oneProject = as.name(paste0("impute_oneProject_", workflow))
  dat_imp_list = lapply(dat_list, impute_oneProject, saveFolder)
  cat("Missing data imputation is finished. \n")
  end_time <- Sys.time()
  print(end_time - start_time)
  return(dat_imp_list)
}

myImputationFun = function(mu, s, minOfProtein, numOfMissing) {
  # used by function impute_one_protein
  output = rtruncnorm(
    n = numOfMissing,
    a = -Inf,
    b = minOfProtein,
    mean = mu,
    sd = s
  )
  return(output)
}


# Pulldown ----------------------------------------------------------------

impute_oneProject_pulldown = function(dat, saveFolder) {
  # For rows with missing raw data:
  #   The imputed value will replace the NA in log_Result_Value column
  #   Repeat and raw_log_Result_Value at this row will remain NA
  this_project = unique(as.character(dat$Project_Experiment_SampleGroup))
  cat(
    "================= Missing data imputation in Project_Experiment_SampleGroup: ",
    this_project,
    " ================= \n"
  )
  
  # one happy gaussian missing data imputation:
  # use all available data to fit a gaussian distribution, and plot
  this_mu = mean(dat$log_Result_Value, na.rm = TRUE)
  this_sd = sd(dat$log_Result_Value, na.rm = TRUE)
  
  # perform imputation for all proteins
  dat_split = split(dat, list(dat$proteinID))
  dat_split_imp = lapply(dat_split, impute_one_protein_pulldown, this_mu, this_sd)
  dat_imp = do.call(rbind, dat_split_imp)
  rm(dat_split, dat_split_imp)
  rownames(dat_imp) <- NULL
  
  # missing proportion table for each sample before and after imputation
  n_protein = length(unique(dat$proteinID))
  cat("\n Missing proportion by SampleName before imputation: ")
  print(round(1 - table(dat$SampleName[!is.na(dat$log_Result_Value)]) /
                n_protein, 3))
  cat("\n Missing proportion by SampleName after imputation: ")
  print(round(1 - table(dat_imp$SampleName[!is.na(dat_imp$log_Result_Value)]) / n_protein, 3))
  return(dat_imp)
}


impute_one_protein_pulldown = function(dat_one_protein, mu, s) {
  # imputation method: one happy gaussian
  missingIDs = is.na(dat_one_protein$log_Result_Value)
  # perform imputation only when missing less than 50% among all wells (SampleNames)
  #   e.g. for one protein, there should be ~96 samples, at least observe 48 data points.
  if (sum(missingIDs) <= 0.5 * nrow(dat_one_protein)) {
    min_i = min(dat_one_protein$log_Result_Value, na.rm = TRUE)
    dat_one_protein$log_Result_Value[missingIDs] =
      myImputationFun(mu, s, minOfProtein = min_i, numOfMissing = sum(missingIDs))
  }
  return(dat_one_protein)
}


# Pro10K ------------------------------------------------------------------

impute_oneProject_pro10K = function(dat, saveFolder) {
  # For rows with missing raw data:
  #   The imputed value will replace the NA in log_Result_Value column
  #   Repeat and raw_log_Result_Value at this row will remain NA
  this_project = unique(as.character(dat$Project_Experiment))
  cat(
    "================= Missing data imputation in Project_Experiment: ",
    this_project,
    " ================= \n"
  )
  if(!("raw_log_Result_Value" %in% names(dat))){
    dat$raw_log_Result_Value = dat$log_Result_Value
  }
  # one happy gaussian missing data imputation:
  # use all available data to fit a gaussian distribution, and plot
  this_mu = mean(dat$log_Result_Value, na.rm = TRUE)
  this_sd = sd(dat$log_Result_Value, na.rm = TRUE)
  
  p1 = ggplot(dat %>% filter(!is.na(log_Result_Value)), aes(log_Result_Value)) +
    geom_histogram(
      aes(y = ..density..),
      binwidth = 0.1,
      fill = "blue",
      alpha = 0.1
    ) +
    geom_density(aes(y = ..density..), size = 1, colour = "blue") +
    stat_function(
      fun = dnorm,
      args = list(
        mean = mean(dat$log_Result_Value, na.rm = TRUE),
        sd = sd(dat$log_Result_Value, na.rm = TRUE)
      ),
      colour = "black",
      size = 1
    ) +
    ggtitle(paste0(
      "Observed data (blue) v.s. Fitted Gaussian distribution (black)"
    ))
  
  # plot observed log_Result_Value v.s. count by protein
  dat_protein = dat %>% filter(!is.na(log_Result_Value)) %>%
    group_by(proteinID) %>%
    summarise(count = n(),
              protein_median = median(log_Result_Value))
  p2 = ggplot(dat_protein, aes(count, protein_median)) +
    geom_point() + geom_smooth(method = "lm") +
    ggtitle("Observed SampleName Counts v.s. median_log_Result_Value for all proteins")
  
  savePlotName = file.path(saveFolder,
                           paste0("Project_", this_project, "_Imputation.jpg"))
  cat(
    paste0(
      "Plot the observed distribution and missing pattern to file: ",
      savePlotName,
      "\n"
    )
  )
  jpeg(filename = savePlotName, width = 1600, height = 800, res = 100, quality = 80)
  gridExtra::grid.arrange(p1, p2, nrow = 1)
  dev.off()
  rm(p1, p2)
  
  # perform imputation for all proteins
  dat_split = split(dat, list(dat$proteinID))
  dat_split_imp = lapply(dat_split, impute_one_protein_pro10K, this_mu, this_sd)
  dat_imp = do.call(rbind, dat_split_imp)
  rm(dat_split, dat_split_imp)
  rownames(dat_imp) <- NULL
  
  # missing proportion table and plot for each well before and after imputation
  n_protein = length(unique(dat$proteinID))
  
  dat_well = dat_imp %>% group_by(SampleName) %>%
    summarise(
      NumProteins_raw = sum(!is.na(raw_log_Result_Value)),
      NumProteins_imp = sum(!is.na(log_Result_Value)),
      outlier = unique(outlier)
    )
  dat_well_long = reshape2::melt(dat_well,
                                 id.var = c('SampleName', 'outlier'),
                                 value.name = "NumProteins")
  colnames(dat_well_long)[which(colnames(dat_well_long) == "variable")] = "Group"
  
  p = ggplot(dat_well_long, aes(SampleName, NumProteins, colour = Group)) +
    geom_point(size = 3, aes(shape = outlier)) +
    geom_hline(
      yintercept = n_protein,
      linetype = "dashed",
      color = "orange",
      size = 1.5
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Number of proteins quantified by well: Before and After Imputation")
  savePlotName = file.path(saveFolder,
                           paste0("Project_", this_project, "_CountProteins_Imp.jpg"))
  cat(
    paste0(
      "Plot the NumProteins by well before and after Imputation to file: ",
      savePlotName,
      "\n"
    )
  )
  jpeg(filename = savePlotName, width = 1200, height = 800, res = 100, quality = 80)
  print(p)
  dev.off()
  rm(p)
  return(dat_imp)
}

impute_one_protein_pro10K = function(dat_one_protein, mu, s) {
  # imputation method: one happy gaussian
  missingIDs = is.na(dat_one_protein$log_Result_Value)
  # perform imputation only when missing less than 50% among all wells (SampleNames)
  #   e.g. for one protein, there should be ~96 samples, at least observe 48 data points.
  if (sum(missingIDs) <= 0.5 * nrow(dat_one_protein)) {
    min_i = min(dat_one_protein$log_Result_Value, na.rm = TRUE)
    dat_one_protein$log_Result_Value[missingIDs] =
      myImputationFun(mu, s, minOfProtein = min_i, numOfMissing = sum(missingIDs))
  }
  return(dat_one_protein)
}


# CETSA -------------------------------------------------------------------

impute_oneProject_CETSA = function(dat, saveFolder) {
  this_project = unique(as.character(dat$Project_Experiment))
  cat(
    "================= Check missing data in Project",
    this_project,
    " ================= \n"
  )
  # Assuming that the missing value only exists in "log_Result_Value" column
  n_missing = sum(is.na(dat$log_Result_Value))
  p_missing = 100 * round(n_missing / nrow(dat), 4)
  if (n_missing > 0) {
    cat(
      paste0(
        "\t Number of missing: " ,
        n_missing ,
        ", \t",
        "Percentage of missing: ",
        p_missing,
        "% \n"
      )
    )
  } else {
    cat("\t No missing record in log_Result_Value.\n")
    return(dat)
  }
  cat(
    "================= Missing data imputation in Project ",
    this_project,
    " ================= \n"
  )
  # one happy gaussian missing data imputation:
  #   for each temperature, use available data to fit a gaussian distribution
  t1 = aggregate(dat$log_Result_Value, by=list(dat$Temperature), FUN=mean, na.rm=TRUE)
  colnames(t1)[which(names(t1) == "Group.1")] = "Temperature"
  colnames(t1)[which(names(t1) == "x")] = "T_avg"
  t2 = aggregate(dat$log_Result_Value, by=list(dat$Temperature), FUN=sd, na.rm=TRUE)
  colnames(t2)[which(names(t2) == "Group.1")] = "Temperature"
  colnames(t2)[which(names(t2) == "x")] = "T_sd"
  tt = merge(t1, t2, by = "Temperature")
  dat = merge(dat, tt, by = "Temperature", all.x = TRUE)
  
  dat_split = split(dat, list(dat$proteinID, dat$Temperature))
  dat_split_imp = lapply(dat_split, impute_one_protein_CETSA)
  dat_imp = do.call(rbind, dat_split_imp)
  rm(dat_split, dat_split_imp)
  rownames(dat_imp) <- NULL
  
  # missing proportion table for each sample at each temperature before and after imputation
  n_protein = length(unique(dat$proteinID))
  cat("\n Missing proportion by SampleGroup and Temperature before imputation: ")
  print(round(1 - table(dat$SampleName[!is.na(dat$log_Result_Value)], 
                        dat$Temperature[!is.na(dat$log_Result_Value)]) / n_protein, 3))
  cat("\n Missing proportion by SampleGroup and Temperature after imputation: ")
  print(round(1 - table(dat_imp$SampleName[!is.na(dat_imp$log_Result_Value)], 
              dat_imp$Temperature[!is.na(dat_imp$log_Result_Value)]) /n_protein, 3))
  return(dat_imp)
}

impute_one_protein_CETSA = function(dat_one_protein) {
  # imputation method: one happy gaussian
  missingIDs = is.na(dat_one_protein$log_Result_Value)
  # perform imputation only when missing less than 50% at each temperature
  #   e.g. for one protein, there should be 6 samples at each temperature, at least observe 3 data points.
  if (sum(missingIDs) <= 0.5 * nrow(dat_one_protein)) {
    min_i = min(dat_one_protein$log_Result_Value, na.rm = TRUE)
    mu = dat_one_protein$T_avg[1]
    s = dat_one_protein$T_sd[1]
    dat_one_protein$log_Result_Value[missingIDs] = 
      myImputationFun(mu, s, minOfProtein = min_i, numOfMissing = sum(missingIDs))
  }
  dat_one_protein$T_avg = NULL
  dat_one_protein$T_sd = NULL
  return(dat_one_protein)
}

impute_oneProject_2dCETSA = impute_oneProject_CETSA
