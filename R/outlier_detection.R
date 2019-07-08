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



## outlier detection

utils::globalVariables(c("."))

#   add a column of outlier flag, but don't remove any data
#   call function flag_outlier, which maybe different for different workflow
flag_outlier_list = function(workflow, dat_list, saveFolder) {
  start_time <- Sys.time()
  # re-define flag_outlier function according to workflow
  #flag_outlier = eval(paste0("flag_outlier_",attr(dat_list,"workflow")))
  #flag_outlier = eval(paste0("flag_outlier_", workflow))
  flag_outlier = as.name(paste0("flag_outlier_", workflow))
  dat_list = lapply(dat_list, flag_outlier, saveFolder)
  cat("Outlier detection (add outlier flag, do not remove data) is finished. \n")
  end_time <- Sys.time()
  print(end_time - start_time)
  return(dat_list)
}

# common function used by some flag_outlier()
small_outlier = function(x, threshold = 5) {
  # Outlier detection method for small sample
  # Z*-score: It differs from Z-scores by replacing mean with median and SD with mad.
  # x is a vector of whatever
  output = abs((x - median(x)) / mad(x))
  output = output > threshold
  return(output)
}

wells_corr = function(dat){
  temp = dat %>% select(proteinID, SampleName, log_Result_Value) %>% 
    spread(key = SampleName, value = log_Result_Value)
  M = cor(temp[,-1], use = "pairwise.complete.obs")
  avgCorWithOthers = (colSums(M, na.rm = TRUE)-1)/(nrow(M)-1)
  output = data.frame(SampleName = names(avgCorWithOthers),avgCorWithOthers)
  return(output)
}

# Pulldown ----------------------------------------------------------------

# flag outlier according to median log_Result_Value OR corr for each SampleName
flag_outlier_pulldown = function(dat, saveFolder){
  this_project = unique(as.character(dat$Project_Experiment_SampleGroup))
  cat("================= Outlier detection for Project_Experiment_SampleGroup: ",this_project," ================= \n")
  stat_SampleName = dat %>% 
    select(proteinID, SampleName, log_Result_Value) %>% 
    filter(!is.na(log_Result_Value)) %>% 
    group_by(SampleName) %>% 
    summarise(median_log_Result_Value = median(log_Result_Value),
              count_proteinID = n()) %>%
    mutate(outlier_numProteins = (small_outlier(median_log_Result_Value) & small_outlier(count_proteinID)))
  #--------------- calculate outlier_corr
  df_temp = wells_corr(dat) %>% mutate(SampleName = as.character(SampleName))
  stat_SampleName = stat_SampleName %>% left_join(df_temp, by = "SampleName") %>% mutate(outlier_corr = avgCorWithOthers < 0.8)
  #--------------- add outlier by OR
  stat_SampleName$outlier = (stat_SampleName$outlier_numProteins | stat_SampleName$outlier_corr)
  #--------------- save outlier statistics to file
  saveFileName = file.path(saveFolder, "outlier_summary.csv")
  stat_SampleName$Project_Experiment_SampleGroup = this_project
  if(file.exists(saveFileName)){
    write.table(stat_SampleName, file = saveFileName, sep = ",", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE, append = TRUE)
  }else{
    write.table(stat_SampleName, file = saveFileName, sep = ",", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
  }
  print(knitr::kable(stat_SampleName))
  stat_SampleName = stat_SampleName %>% select(SampleName, outlier_numProteins, outlier_corr, outlier)
  dat = dat %>% left_join(stat_SampleName, by = "SampleName")
  return(dat)
}

# Pro10K ------------------------------------------------------------------
# Outlier Detection at well(SampleName) level
# find outlier for data from one project_expriment
# Add three columns to the dataframe:
#   outlier_numProteins
#   outlier_corr
#   outlier = outlier_numProteins | outlier_corr

flag_outlier_pro10K = function(dat, saveFolder){
  this_project = unique(as.character(dat$Project_Experiment))
  cat("================= Outlier detection for Project_Experiment: ",this_project," ================= \n")
  
  dat_well = as_tibble(dat) %>% filter(!is.na(log_Result_Value)) %>% 
    group_by(Project_Experiment, SampleGroup, SampleName) %>% 
    summarise(NumProteins = n())
  
  dat_plate = dat_well %>%
    group_by(Project_Experiment) %>% 
    summarise(Q1 = quantile(NumProteins, 0.25), 
              Q3 = quantile(NumProteins, 0.75)) %>%
    mutate(lowerCutoff = Q1-1.5*(Q3-Q1)) 
  
  dat_well = dat_well %>%
    left_join(dat_plate, by = "Project_Experiment") %>%
    mutate(outlier_numProteins = NumProteins < lowerCutoff) %>%
    select(-c(Q1, Q3,lowerCutoff)) 
  # calculate outlier_corr
  df_temp = wells_corr(dat) %>% mutate(SampleName = as.character(SampleName))
  dat_well = dat_well %>% left_join(df_temp, by = "SampleName") %>% mutate(outlier_corr = avgCorWithOthers < 0.8)
  # add outlier by OR
  dat_well$outlier = (dat_well$outlier_numProteins | dat_well$outlier_corr)
  #--------------- save outlier statistics to file
  saveFileName = file.path(saveFolder, "outlier_summary.csv")
  dat_well$Project_Experiment = this_project
  if(file.exists(saveFileName)){
    write.table(dat_well, file = saveFileName, sep = ",", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE, append = TRUE)
  }else{
    write.table(dat_well, file = saveFileName, sep = ",", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
  }
  # print outlier wells (if any):
  if(sum(dat_well$outlier) > 0){
    cat("Outliers: \n")
    print(knitr::kable(dat_well %>% filter(outlier)))
  }else{
    cat("No outliers found. \n")
  }
  #--------------- save outlier plot
  p = ggplot(dat_well, aes(SampleName, NumProteins)) +
    geom_point(size =3, aes(shape = outlier_numProteins, colour = outlier_corr)) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    ggtitle("Number of proteins quantified by well")
  savePlotName = file.path(saveFolder,
                           paste0("Project_",this_project,"_Outliers.jpg")
  )
  cat(paste0("Plot the NumProteins by well to file: ",savePlotName,"\n"))
  jpeg(filename = savePlotName,
       width = 1200, height = 800,
       res = 100,
       quality = 80)
  print(p)
  dev.off()
  rm(p)
  #--------------- add the three columns of outlier flag (TRUE/FALSE) to the original dat 
  dat = as_tibble(dat) %>% 
    left_join(dat_well %>% select(-NumProteins, -avgCorWithOthers), 
              by = c("Project_Experiment","SampleGroup","SampleName"))
  return(dat)
}


# CETSA -------------------------------------------------------------------

# flag outlier according to median log_Result_Value for each SampleName and wellCorr
flag_outlier_CETSA = function(dat, saveFolder){
  this_project = unique(as.character(dat$Project_Experiment))
  cat("================= Outlier detection for Project_Experiment_SampleGroup: ",this_project," ================= \n")
  stat_SampleName = dat %>% 
    select(proteinID, Temperature, SampleName, log_Result_Value) %>% 
    filter(!is.na(log_Result_Value)) %>% 
    group_by(Temperature, SampleName) %>% 
    summarise(median_log_Result_Value = median(log_Result_Value),
              count_proteinID = n()) %>%
    mutate(outlier_numProteins = (small_outlier(median_log_Result_Value) | small_outlier(count_proteinID)))
  # calculate outlier_corr
  temp_list = dat %>% select(proteinID, Temperature, SampleName, log_Result_Value)
  temp_list = split(temp_list, temp_list$Temperature)
  temp_list_2 = lapply(temp_list, wells_corr) # add avgCorWithOthers column
  for(i in 1:length(temp_list_2)){
    temp_list_2[[i]]$Temperature = names(temp_list_2)[i]
  }
  temp = do.call(rbind, temp_list_2)
  temp$Temperature = as.numeric(temp$Temperature)
  temp$SampleName = as.character(temp$SampleName)
  stat_SampleName = stat_SampleName %>% 
    inner_join(temp, by = c("Temperature", "SampleName")) %>% 
    mutate(outlier_corr = avgCorWithOthers < 0.8)
  # add outlier by OR
  stat_SampleName$outlier = (stat_SampleName$outlier_numProteins | stat_SampleName$outlier_corr)
  #--------------- save outlier statistics to file
  saveFileName = file.path(saveFolder, "outlier_summary.csv")
  stat_SampleName$Project_Experiment = this_project
  if(file.exists(saveFileName)){
    write.table(stat_SampleName, file = saveFileName, sep = ",", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE, append = TRUE)
  }else{
    write.table(stat_SampleName, file = saveFileName, sep = ",", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
  }
  # plot
  print(knitr::kable(stat_SampleName))
  p1 = ggplot(stat_SampleName, aes(Temperature, median_log_Result_Value, colour = SampleName)) + 
    geom_point(size=3,aes(shape = outlier)) 
  p2 = ggplot(stat_SampleName, aes(Temperature, count_proteinID, colour = SampleName)) + 
    geom_point(size=3,aes(shape = outlier))+ geom_line()
  savePlotName = file.path(saveFolder,paste0("Project_",this_project,"_Outlier.jpg"))
  cat(paste0("Plot the Outlier detection results to file: ",savePlotName,"\n"))
  jpeg(filename = savePlotName,
       width = 1000, height = 1200,
       res = 100,
       quality = 80)
  gridExtra::grid.arrange(p1, p2, ncol=1, nrow = 2)
  dev.off()
  rm(p1,p2)
  
  stat_SampleName = stat_SampleName %>% 
    select(SampleName, Temperature, outlier_numProteins, outlier_corr, outlier)
  dat = dat %>% left_join(stat_SampleName, by = c("SampleName","Temperature"))
  return(dat)
}

flag_outlier_2dCETSA = flag_outlier_CETSA
