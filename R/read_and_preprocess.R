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



## Read raw data and preprocess:
## - log transform intensity data
## - add some columns needed for analysis
## - expand data so all samplegroup by protein combination are present (if missing)

# Pulldown ----------------------------------------------------------------

readRawData_pulldown = function(dataPath,
                                resultType = "Intensity",
                                baseline = "DMSO") {
  cat("================= Read and format raw data ================= \n")
  start_time <- Sys.time()
  
  dat_raw = read.table(
    dataPath,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  # choose subset of rows where Result_Type is the resultType value
  #   change all to lower case
  resultType = tolower(resultType)
  dat_raw$Result_Type = tolower(dat_raw$Result_Type)
  if (! resultType %in% unique(dat_raw$Result_Type)) {
    stop(paste0("Result_Type ", resultType, " not found in the input data"))
  }
  dat_raw = dat_raw[dat_raw$Result_Type == resultType, ]
  
  # check whether the user input baseline exists in Dose
  if (! baseline %in% (dat_raw$Dose)) {
    stop(paste0("baseline ", baseline, " not found in the input data"))
  }
  
  # combine ProjectName, ExperimentName and SampleGroup as a factor to split data
  dat_raw$Project_Experiment_SampleGroup = paste(dat_raw$ProjectName,
                                                 dat_raw$ExperimentName,
                                                 dat_raw$SampleGroup,
                                                 sep = "+")
  
  # Keep useful columns in raw data
  #dat = dat_raw[, c("proteinID", "SampleName", "Repeat",
  #                  "Dose", "Result_Value","Project_Experiment_SampleGroup")]
  dat = dat_raw
  
  # split the dataframe into a list of dataframes according to Project_Experiment_SampleGroup column
  # if there is only one factor, it will be a list with only one element
  all_projectExpr = unique(dat$Project_Experiment_SampleGroup)
  cat(
    length(all_projectExpr),
    "projects & experiments & samplegroups: ",
    all_projectExpr,
    "\n"
  )
  cat("Split the raw data according to Project_Experiment_SampleGroup factor.\n")
  dat_list = split(dat, f = dat$Project_Experiment_SampleGroup)
  
  # Pre-process/format the datasets by Project_Experiment
  dat_list = lapply(dat_list, preprocessRawData_pulldown, baseline)
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  return(dat_list)
}

preprocessRawData_pulldown = function(dat, baseline) {
  # This function is used by function readRawData.
  this_project = unique(dat$Project_Experiment_SampleGroup)
  cat(
    "\n ------------ Preprocess data for Project_Experiment_SampleGroup:",
    this_project,
    "------------ \n"
  )
  # Take log of the Result_Value
  #cat("Take log 10 of Result_Value for saving to dat_processed.csv \n")
  #dat$log10_Result_Value = log10(dat$Result_Value)
  cat("Take log 2 of Result_Value for model fitting. \n")
  dat$log_Result_Value = log2(dat$Result_Value)
  
  # Change data type from chr to factor for proteinID
  dat$proteinID = as.factor(dat$proteinID)
  # order the Condition factor according to its corresponding Dose
  #temp = as_tibble(unique(dat[,c("Dose","Condition")])) %>% arrange(Dose)
  dat$Dose_F = factor(dat$Dose, levels = sort(unique(dat$Dose)))
  # make sure baseline of Dose is as specified
  dat$Dose_F <- relevel(dat$Dose_F, as.character(baseline))
  
  # print several summary statistics
  cat("Number of proteins: ", length(unique(dat$proteinID)), "\n")
  all_Dose = sort(unique(dat$Dose))
  cat(length(all_Dose), "Doses: ", all_Dose, "\n")
  cat("Number of SampleNames: ", length(unique(dat$SampleName)), "\n")
  
  # Extend the data frame so that missing data will be easily identify
  # complete data should have all combinations of proteinID and SampleName
  allgroup_df = merge(unique(dat["proteinID"]),
                      unique(dat[, c("SampleName",
                                     "Dose",
                                     "Dose_F",
                                     "Repeat",
                                     "Project_Experiment_SampleGroup")]),
                      all = TRUE)
  #stopifnot(nrow(allgroup_df) == 
  #            length(unique(dat$proteinID)) * length(unique(dat$SampleName)))
  dat = merge(dat,
              allgroup_df,
              by = names(allgroup_df),
              all.y = TRUE)
  stopifnot(nrow(allgroup_df) == nrow(dat))
  
  # Now the missing value only exists in "log_Result_Value" column
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
  } else{
    cat("\t No missing record in log_Result_Value. \n ")
  }
  
  # add a column of numeric logDose for model fitting
  dat = add_logDose(dat)
  
  return(dat)
}



# Pro10K ------------------------------------------------------------------

readRawData_pro10K = function(dataPath,
                              resultType = "Intensity",
                              baseline = "DMSO") {
  cat("================= Read and format raw data ================= \n")
  start_time <- Sys.time()
  
  cat(paste0("Read raw data from file: ", dataPath, "\n"))
  
  dat_raw = read.table(
    dataPath,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  # choose subset of rows where Result_Type is the resultType value
  dat_raw = dat_raw[dat_raw$Result_Type == resultType,]
  
  # check whether the user input baseline exists in SampleGroup levels
  if (!baseline %in% (dat_raw$SampleGroup)) {
    stop(paste0(
      "The SampleGroup baseline ",
      baseline,
      " doesn't exists in the input data. "
    ))
  }
  
  # combine ProjectName and ExperimentName as a factor to split data
  dat_raw$Project_Experiment = paste(dat_raw$ProjectName, dat_raw$ExperimentName, sep = "+")
  
  # Keep useful columns in raw data
  #dat = dat_raw[, c("proteinID", "SampleGroup","SampleName", "Repeat",
  #                  "Result_Value","Project_Experiment")]
  dat = dat_raw
  
  # split the dataframe into a list of dataframes according to ProjectName column
  # if there is only one ProjectName, it will be a list with only one element
  all_projectExpr = unique(dat$Project_Experiment)
  cat(length(all_projectExpr),
      "projects & experiments: ",
      all_projectExpr,
      "\n")
  cat("Split the raw data according to Project_Experiment factor.\n")
  dat_list = split(dat, f = dat$Project_Experiment)
  
  # Pre-process/format the datasets by Project_Experiment
  dat_list = lapply(dat_list, preprocessRawData_pro10K, baseline)
  
  end_time <- Sys.time()
  print(end_time - start_time)
  return(dat_list)
}

preprocessRawData_pro10K = function(dat, baseline) {
  # This function is used by function readRawData.
  this_project = unique(dat$Project_Experiment)
  cat(
    "\n ------------ Preprocess data for Project_Experiment:",
    this_project,
    "------------ \n"
  )
  # Check whether there is any missing
  #if(!all(complete.cases(dat))){
  #  print(dat[!complete.cases(dat),])
  #  stop("Unexpected missing cases in the input data. ")
  #}
  
  # Take log of the Result_Value
  #cat("Take log 10 of Result_Value for saving to dat_processed.csv \n")
  #dat$log10_Result_Value = log10(dat$Result_Value)
  cat("Take log 2 of Result_Value for model fitting. \n")
  dat$log_Result_Value = log2(dat$Result_Value)
  
  # Change data type from chr to factor for SampleGroup, SampleName, proteinID
  dat$SampleGroup = as.factor(dat$SampleGroup)
  dat$SampleGroup <- relevel(dat$SampleGroup, baseline)
  dat$proteinID = as.factor(dat$proteinID)
  
  # print several summary statistics
  cat("Number of proteins: ", length(unique(dat$proteinID)), "\n")
  all_SampleGroup = levels(dat$SampleGroup)
  cat(length(all_SampleGroup),
      "Sample Groups: ",
      all_SampleGroup,
      "\n")
  cat("Number of SampleNames (Wells): ", length(unique(dat$SampleName)), "\n")
  sampleGroup_repeat = as_tibble(dat) %>% filter(!(SampleGroup %in% c(baseline, "Empty"))) %>%
    group_by(SampleGroup) %>% summarise(numberOfRepeat = max(Repeat))
  cat("Table of Repeats for non-baseline & non-empty SampleGroup: ")
  print(knitr::kable(sampleGroup_repeat))
  
  # Extend the data frame so that missing data will be easily identify
  # it needs to be done after split by Project/Expr since the SampleGroup is different.
  # complete data should have all combinations of proteinID and SampleName
  allgroup_df = merge(unique(dat["proteinID"]),
                      unique(dat[, c("SampleName", "SampleGroup", "Project_Experiment", "Repeat")]),
                      all = TRUE)
  #stopifnot(nrow(allgroup_df) == length(unique(dat$proteinID)) * length(unique(dat$SampleName)))
  dat = merge(
    dat,
    allgroup_df,
    by = c(
      "proteinID",
      "SampleGroup",
      "SampleName",
      "Project_Experiment",
      "Repeat"
    ),
    all.y = TRUE
  )
  stopifnot(nrow(allgroup_df) == nrow(dat))
  
  # Now the missing value only exists in "log_Result_Value" column
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
  } else{
    cat("\t No missing record in log_Result_Value. \n ")
  }
  
  return(dat)
}

# CETSA -------------------------------------------------------------------

readRawData_CETSA = function(dataPath,
                             resultType = "Intensity",
                             baseline = "DMSO") {
  cat("================= Read and format raw data ================= \n")
  start_time <- Sys.time()
  
  cat(paste0("Read raw data from file: ", dataPath, "\n"))
  dat_raw = read.table(
    dataPath,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  # choose subset of rows where Result_Type is the resultType value
  if (!resultType %in% unique(dat_raw$Result_Type)) {
    stop(paste0(
      "Result_Type ",
      resultType,
      " doesn't exists in the input data. "
    ))
  }
  dat_raw = dat_raw[dat_raw$Result_Type == resultType,]
  
  # check whether the user input baseline exists in SampleGroup levels
  if (!baseline %in% (dat_raw$SampleGroup)) {
    stop(paste0(
      "The SampleGroup baseline ",
      baseline,
      " doesn't exists in the input data. "
    ))
  }
  #dat_raw$SampleGroup = as.factor(dat_raw$SampleGroup)
  #stopifnot(baseline %in% levels(dat_raw$SampleGroup))
  #dat_raw$SampleGroup<- relevel(dat_raw$SampleGroup, baseline)
  
  # combine ProjectName and ExperimentName as a factor to split data
  dat_raw$Project_Experiment = paste(dat_raw$ProjectName, dat_raw$ExperimentName, sep = "+")
  
  # Keep useful columns in raw data
  #dat = dat_raw[, c("proteinID", "SampleGroup","SampleName", "Repeat","Temperature","Dose",
  #                  "Result_Value","Project_Experiment")]
  dat = dat_raw
  
  # split the dataframe into a list of dataframes according to Project_Experiment column
  # if there is only one Project_Experiment, it will be a list with only one element
  all_projectExpr = unique(dat$Project_Experiment)
  cat(length(all_projectExpr),
      "projects & experiments: ",
      all_projectExpr,
      "\n")
  cat("Split the raw data according to Project_Experiment factor.\n")
  dat_list = split(dat, f = dat$Project_Experiment)
  
  # Pre-process/format the datasets by Project_Experiment
  dat_list = lapply(dat_list, preprocessRawData_CETSA, baseline)
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  return(dat_list)
}

preprocessRawData_CETSA = function(dat, baseline) {
  # This function is used by function readRawData.
  this_project = unique(dat$Project_Experiment)
  cat(
    "\n ------------ Preprocess data for Project_Experiment:",
    this_project,
    "------------ \n"
  )
  #dat = dat[complete.cases(dat),] #the useless columns may be NAs
  
  # Take log of the Result_Value
  #cat("Take log 10 of Result_Value for saving to dat_processed.csv \n")
  #dat$log10_Result_Value = log10(dat$Result_Value)
  cat("Take log 2 of Result_Value for model fitting. \n")
  dat$log_Result_Value = log2(dat$Result_Value)
  
  # Change data type from chr to factor for SampleGroup, proteinID
  dat$SampleGroup = factor(dat$SampleGroup, levels = unique(dat[order(dat$Dose), "SampleGroup"]))
  dat$SampleGroup <- relevel(dat$SampleGroup, baseline)
  dat$proteinID = as.factor(dat$proteinID)
  
  # print several summary statistics
  all_temperatures = sort(unique(dat$Temperature))
  all_samples = sort(as.character(unique(dat$SampleName)))
  all_proteins = sort(as.character(unique(dat$proteinID)))
  all_samplegroups = sort(as.character(unique(dat$SampleGroup)))
  cat("Number of proteins: ", length(all_proteins), "\n")
  cat(length(all_samplegroups),
      "SampleGroups: ",
      all_samplegroups,
      "\n")
  cat(length(all_samples), "SampleNames: ", all_samples, "\n")
  cat(length(all_temperatures),
      "temperatures: ",
      all_temperatures,
      "\n")
  
  # Extend the data frame so that missing data will be easily identify
  # it needs to be done after split by Project/Expr since the SampleGroup is different.
  # complete data should have all combinations of proteinID, SampleName, Temperature
  allgroup_df = expand.grid(all_proteins, all_temperatures, all_samples)
  colnames(allgroup_df) = c("proteinID", "Temperature", "SampleName")
  allgroup_df_2 = unique(dat[, c("SampleGroup", "Dose", "SampleName", "Project_Experiment")])
  allgroup_df = merge(allgroup_df, allgroup_df_2, by = "SampleName", all.x = TRUE)
  dat = merge(
    dat,
    allgroup_df,
    by = c(
      "proteinID",
      "Temperature",
      "SampleGroup",
      "Dose",
      "SampleName",
      "Project_Experiment"
    ),
    all.y = TRUE
  )
  stopifnot(nrow(allgroup_df) == nrow(dat))
  
  return(dat)
}

readRawData_2dCETSA = readRawData_CETSA
