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



## main function of the rCPDMS package

MSanalysis = function(workflow,
                  dataPath,
                  resultType = "Intensity",
                  baseline = "DMSO",
                  saveFolder = "output",
                  outlier = FALSE,
                  normalization = TRUE,
                  imputation = FALSE) {
  # check whether workflow is valid
  if (!(workflow %in% c("pulldown", "pro10K", "CETSA", "2dCETSA"))) {
    stop(paste0("Input parameter workflow = ", workflow,
                " is not available."))
  }

  # check whether dataPath is valid
  if (!file.exists(dataPath)) {
    stop(paste0("dataPath ", dataPath, " not found"))
  }
  
  # check whether the specified saveFolder exists, otherwise create new folder
  if (!file.exists(saveFolder)) {
    if (! dir.create(saveFolder, recursive = TRUE)) {
        stop(paste0("error creating the folder ", saveFolder))
    }
  }
  
  # Read and format raw data
  dat_list = do.call(paste0("readRawData_", workflow),
                     args = list(dataPath, resultType, baseline))
  
  # Outlier detection according to missing number of proteins
  if (outlier) {
    dat_list = flag_outlier_list(workflow, dat_list, saveFolder)
  } else {
    cat("------- No outlier detection -------\n")
    for(dat in dat_list){dat$outlier = FALSE}
  }
  
  # Normalization
  if (normalization) {
    dat_list = normalize_list(workflow, dat_list, saveFolder)
  } else {
    cat("------- No normalization -------\n")
    for(dat in dat_list){dat$raw_log_Result_Value = dat$log_Result_Value}
  }
  
  # Missing data imputation
  if (imputation) {
    dat_list = impute_list(workflow, dat_list, saveFolder)
  } else {
    cat("------- No missing data imputation -------\n")
  }
  
  # Save pre-processed data (combined across Projects/Experiments/SampleGroups)
  saveProcessedData_list(dat_list, saveFolder)
  
  # Fitting model for all projects
  results_list = do.call(paste0("fitAll_", workflow), args = list(dat_list, baseline))
  
  # Add FDR correction to p-values in results
  results_list = do.call(paste0("add_FDR_", workflow), args = list(results_list, baseline))
  
  # Finally, combine results from multiple projects/experiments and save one results dataframe to file:
  results = saveResults_list(results_list, saveFolder)
  
  # Part 2 of analysis for 2dCETSA: Do the same analysis as pulldown for each temperature
  if (workflow == "2dCETSA") {
    results_list2 = fitAll_2dCETSA_by_temp(dat_list, baseline)
    results_list2 = add_FDR_2dCETSA_by_temp(results_list2)
    # Save results_list2 to file:
    results2 = saveResults_list(results_list2, saveFolder, fileName = "results2")
  }
  return(results)
}