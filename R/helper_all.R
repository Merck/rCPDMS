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



## misc. utility functions

utils::globalVariables(c("Dose"))

# Used in pulldown and 2dCETSA
add_logDose = function(dat) {
  #df_Dose = dat %>% distinct(Dose, Dose_F) %>% mutate(logDose = log10(Dose)) 
  df_Dose = dat %>% distinct(Dose) %>% mutate(logDose = log10(Dose)) 
  
  ## Make DMSO one log step below the lowest concentration
  logPosDose <- log10(sort(unique(dat$Dose[dat$Dose > 0])))
  if (length(logPosDose)>1) {
    df_Dose$logDose[df_Dose$Dose == 0] <- logPosDose[1] - (logPosDose[2] - logPosDose[1])
  }
  #dat = dat %>% left_join(df_Dose, by = c("Dose","Dose_F"))
  dat = dat %>% left_join(df_Dose, by = "Dose")
  return(dat)
}

# save processed dataframe to file
saveProcessedData_list = function(dat_list, saveFolder){
  dat_toSave = do.call(rbind, dat_list)
  cat("Calculate log10_Result_Value for saving to dat_processed.csv \n")
  dat_toSave$log10_Result_Value = sapply(dat_toSave$log_Result_Value, log2_to_log10)
  if("raw_log_Result_Value" %in% colnames(dat_toSave)){
    cat("Calculate raw_log10_Result_Value for saving to dat_processed.csv \n")
    dat_toSave$raw_log10_Result_Value = sapply(dat_toSave$raw_log_Result_Value, log2_to_log10)
  }
  saveFileName = file.path(saveFolder,"dat_processed.csv")
  cat(paste0("Save pre-processed dataset to one CSV file: ",saveFileName,"\n"))
  write.csv(dat_toSave, file = saveFileName, row.names = F, na = "")
}

log2_to_log10 = function(x){
  return(x*log10(2))
}

# save results dataframe to file
saveResults_list = function(results_list, saveFolder, fileName = "results"){
  results = bind_rows(results_list)
  if("Project_Experiment_SampleGroup" %in% colnames(results)){
    results = results %>% arrange(Project_Experiment_SampleGroup, proteinID)
  }else{
    if("SampleGroup" %in% colnames(results)){
      results = results %>% arrange(Project_Experiment, proteinID, SampleGroup)
    }else{
      results = results %>% arrange(Project_Experiment, proteinID)
    }
  }
  saveFileName = file.path(saveFolder, paste0(fileName,".csv"))
  cat(paste0("Save results dataframe to one CSV file: ",saveFileName,"\n"))
  write.csv(results, file = saveFileName, row.names = FALSE, na = "")
  return(results)
}

# erfc is used in significance A, B
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)

## Significance A,B MaxQuant
significance_A = function(r){
  # input r is a vector of log ratios
  r_1 = quantile(r, probs = 0.1587, na.rm = TRUE)
  r0 = quantile(r, probs = 0.5, na.rm = TRUE)
  r1 = quantile(r, probs = 0.8413, na.rm = TRUE)
  z = numeric(length = length(r))
  for(i in 1:length(r)){
    if(!is.na(r[i])){
      if(r[i] >= r0){
        z[i] = (r[i]-r0)/(r1-r0)
      }else{
        z[i] = (r0-r[i])/(r0-r_1)
      }
    }else{
      z[i] = NA
    }
  }
  output = as.numeric(0.5*erfc(z/sqrt(2)))
  return(output)
}

significance_B = function(df){
  r_ind = grep("est_logRatioDose_", names(df))
  new_col = gsub("est_","pVal_sigB_",names(df)[r_ind])
  stopifnot("group" %in% names(df))
  df_split = split(df, df$group)
  for(i in 1:length(df_split)){
    df_one = df_split[[i]]
    df_one[new_col] = significance_A(unlist(df_one[,r_ind]))
    df_split[[i]] = df_one
  }
  output = do.call(rbind, df_split)
  return(output)
}
