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



## Data normalization

# call function normalize_oneProject, which maybe different for different workflow
normalize_list = function(workflow, dat_list, saveFolder = "output") {
  op <- options(warn=-1)
  on.exit(options(op))
  start_time <- Sys.time()
  # re-define normalize_oneProject function according to workflow
  #normalize_oneProject = eval(paste0("normalize_oneProject_", workflow))
  normalize_oneProject = as.name(paste0("normalize_oneProject_", workflow))
  print(normalize_oneProject)
  dat_norm_list = lapply(dat_list, normalize_oneProject, saveFolder)
  cat("Normalization is finished. \n")
  end_time <- Sys.time()
  print(end_time - start_time)
  return(dat_norm_list)
}

# Pulldown ----------------------------------------------------------------

# match the median across SampleName
# normalization for data from one project_expriment_samplegroup
normalize_oneProject_pulldown = function(dat, saveFolder) {
  this_project = unique(as.character(dat$Project_Experiment_SampleGroup))
  cat(
    "================= Normalize data in Project_Experiment_SampleGroup: ",
    this_project,
    " ================= \n"
  )
  # Step 1: Normalize across good wells (SampleName) 
  if(!("outlier" %in% names(dat))){
    # if didn't perform outlier detection before, then treated all wells as good wells.
    dat$outlier = FALSE 
  }
  #   Summary statistics by SampleName
  dat_well = as_tibble(dat) %>% filter(outlier == FALSE) %>%
    group_by(Project_Experiment_SampleGroup, Dose_F, SampleName) %>%
    summarise(median_log_Result_Value = median(log_Result_Value, na.rm = TRUE))
  
  dat_plate = dat_well %>%
    group_by(Project_Experiment_SampleGroup) %>%
    summarise(plate_median = median(median_log_Result_Value))
  
  dat_well = dat_well %>%
    left_join(dat_plate, by = "Project_Experiment_SampleGroup") %>%
    mutate(delta = median_log_Result_Value - plate_median)
  
  dat_norm = as_tibble(dat) %>%
    left_join(
      dat_well %>% dplyr::select(Project_Experiment_SampleGroup, Dose_F, SampleName, delta),
      by = c("Dose_F", "SampleName", "Project_Experiment_SampleGroup")
    ) %>%
    mutate(raw_log_Result_Value = log_Result_Value) %>%
    mutate(log_Result_Value = log_Result_Value - delta) %>%
    dplyr::select(-delta)
  
  dat_norm$log_Result_Value[dat_norm$outlier] = dat_norm$raw_log_Result_Value[dat_norm$outlier]
  
  # Step 2: for the bad wells (if any), match the median of subset of common proteins
  if(any(dat$outlier)){
    protein_subsets = dat %>% filter(outlier) %>% filter(!is.na(log_Result_Value)) %>%
      group_by(proteinID) %>% summarise(count = length(unique(SampleName))) %>% 
      filter(count == max(count)) %>% ungroup() %>% select(proteinID)
    
    dat_well_2 = dat_norm %>% 
      filter(proteinID %in% protein_subsets$proteinID) %>% 
      group_by(Project_Experiment_SampleGroup, SampleName) %>% 
      summarise(median_log_Result_Value = median(log_Result_Value, na.rm = TRUE), outlier = unique(outlier))
    
    dat_plate_2 = dat_well_2 %>% filter(outlier==FALSE) %>% 
      group_by(Project_Experiment_SampleGroup) %>% 
      summarise(plate_median = median(median_log_Result_Value))
    
    dat_well_2 = dat_well_2 %>% filter(outlier) %>%
      left_join(dat_plate_2, by = "Project_Experiment_SampleGroup") %>%
      mutate(delta = median_log_Result_Value - plate_median)
    
    dat_norm_2 = as_tibble(dat) %>% filter(outlier) %>%
      left_join(dat_well_2 %>% select(Project_Experiment_SampleGroup, SampleName, delta), 
                by = c("SampleName", "Project_Experiment_SampleGroup")) %>%
      mutate(raw_log_Result_Value = log_Result_Value) %>% 
      mutate(log_Result_Value = log_Result_Value - delta) %>%
      select(-delta)
    
    dat_norm = dat_norm %>% filter(!outlier) %>% rbind(dat_norm_2) %>% arrange(SampleName)
  }
  p1 = ggplot(dat_norm,
              aes(SampleName, raw_log_Result_Value, colour = Dose_F)) +
    geom_boxplot() + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Boxplot of log(Result_Value) by SampleName: Before normalization")
  p2 = ggplot(dat_norm, aes(SampleName, log_Result_Value, colour = Dose_F)) +
    geom_boxplot() + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Boxplot of log(Result_Value) by SampleName: After normalization")
  savePlotName = file.path(saveFolder,
                           paste0("Project_", this_project, "_Normalization.jpg"))
  cat(paste0("Plot the Normalization results to file: ", savePlotName, "\n"))
  jpeg(
    filename = savePlotName,
    width = 1200,
    height = 800,
    res = 100,
    quality = 80
  )
  gridExtra::grid.arrange(p1, p2, ncol = 2, nrow = 1)
  dev.off()
  rm(p1, p2)
  return(dat_norm)
}


# Pro10K ------------------------------------------------------------------

# Normalization for data from one project_expriment
normalize_oneProject_pro10K = function(dat, saveFolder){
  this_project = unique(as.character(dat$Project_Experiment))
  cat("================= Normalize data in Project_Experiment: ",this_project," ================= \n")
  
  # Step 1: Normalize across good wells (SampleName) 
  if(!("outlier" %in% names(dat))){
    # if didn't perform outlier detection before, then treated all wells as good wells.
    dat$outlier = FALSE 
  }
  #   Summary statistics by well (or SampleName)
  dat_well = as_tibble(dat) %>% filter(outlier == FALSE) %>% 
    group_by(Project_Experiment, SampleGroup, SampleName) %>% 
    summarise(median_log_Result_Value = median(log_Result_Value, na.rm = TRUE))
  
  dat_plate = dat_well %>%
    group_by(Project_Experiment) %>% 
    summarise(plate_median = median(median_log_Result_Value))
  
  dat_well = dat_well %>%
    left_join(dat_plate, by = "Project_Experiment") %>%
    mutate(delta = median_log_Result_Value - plate_median)
  
  dat_norm = as_tibble(dat) %>% 
    left_join(dat_well %>% select(Project_Experiment, SampleGroup, SampleName, delta), 
              by = c("SampleGroup", "SampleName", "Project_Experiment")) %>%
    mutate(raw_log_Result_Value = log_Result_Value) %>% 
    mutate(log_Result_Value = log_Result_Value - delta) %>%
    select(-delta)
  dat_norm$log_Result_Value[dat_norm$outlier] = dat_norm$raw_log_Result_Value[dat_norm$outlier]
  
  # Step 2: for the bad wells (if any), match the median of subset of common proteins
  if(any(dat$outlier)){
  protein_subsets = dat %>% filter(outlier) %>% filter(!is.na(log_Result_Value)) %>%
    group_by(proteinID) %>% summarise(count = length(unique(SampleName))) %>% 
    filter(count == max(count)) %>% ungroup() %>% select(proteinID)
  
  dat_well_2 = dat_norm %>% 
    filter(proteinID %in% protein_subsets$proteinID) %>% 
    group_by(Project_Experiment, SampleGroup, SampleName) %>% 
    summarise(median_log_Result_Value = median(log_Result_Value, na.rm = TRUE), outlier = unique(outlier))
  
  dat_plate_2 = dat_well_2 %>% filter(outlier==FALSE) %>% 
    group_by(Project_Experiment) %>% 
    summarise(plate_median = median(median_log_Result_Value))
  
  dat_well_2 = dat_well_2 %>% filter(outlier) %>%
    left_join(dat_plate_2, by = "Project_Experiment") %>%
    mutate(delta = median_log_Result_Value - plate_median)
  
  dat_norm_2 = as_tibble(dat) %>% filter(outlier) %>%
    left_join(dat_well_2 %>% select(Project_Experiment, SampleGroup, SampleName, delta), 
              by = c("SampleGroup", "SampleName", "Project_Experiment")) %>%
    mutate(raw_log_Result_Value = log_Result_Value) %>% 
    mutate(log_Result_Value = log_Result_Value - delta) %>%
    select(-delta)
  
  dat_norm = dat_norm %>% filter(!outlier) %>% rbind(dat_norm_2) %>% arrange(SampleName)
  }
  p1 = ggplot(dat_norm, aes(SampleName, raw_log_Result_Value, colour = SampleGroup)) +
    geom_boxplot() + theme(legend.position="none") +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    ggtitle("Boxplot of log(Result_Value) by SampleName: Before normalization")
  p2 = ggplot(dat_norm, aes(SampleName, log_Result_Value, colour = SampleGroup)) +
    geom_boxplot() + theme(legend.position="none") +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    ggtitle("Boxplot of log(Result_Value) by SampleName: After normalization")
  savePlotName = file.path(saveFolder,
                           paste0("Project_",this_project,"_Normalization.jpg")
  )
  cat(paste0("Plot the Normalization results to file: ",savePlotName,"\n"))
  jpeg(filename = savePlotName,
       width = 1000, height = 1200,
       res = 100,
       quality = 80)
  gridExtra::grid.arrange(p1, p2, ncol=1, nrow = 2)
  dev.off()
  rm(p1,p2)
  return(dat_norm)
}


# CETSA -------------------------------------------------------------------

# normalization for data from one project
normalize_oneProject_CETSA = function(dat, saveFolder = "output") {
  this_project = unique(as.character(dat$Project_Experiment))
  cat("================= Normalize data in Project",
      this_project,
      " ================= \n")
  
  # boxplot of samples/temperature before normalization
  p1 = ggplot(dat[!is.na(dat$log_Result_Value), ], aes(as.factor(Temperature), log_Result_Value)) +
    geom_boxplot(aes(fill = SampleName)) +
    ggtitle("Boxplot by Sample and Temperature - Before Normalization") +
    xlab("Temperature")
  
  # one-step normalization:
  #   Normalize across samples at each temperature
  dat_split_by_Temperature = split(dat, f = dat$Temperature)
  dat = do.call(rbind,
                lapply(dat_split_by_Temperature, normalize_dat_one_T))
  
  # boxplot of samples/temperature after normalization
  p2 = ggplot(dat[!is.na(dat$log_Result_Value), ], aes(as.factor(Temperature), log_Result_Value)) +
    geom_boxplot(aes(fill = SampleName)) +
    ggtitle("Boxplot by Sample and Temperature - After Normalization") +
    xlab("Temperature")
  savePlotName = file.path(saveFolder,
                           paste0("Project_", this_project, "_Normalization.jpg"))
  cat(paste0("Plot the Normalization results to file: ", savePlotName, "\n"))
  jpeg(
    filename = savePlotName,
    width = 1000,
    height = 1200,
    res = 100,
    quality = 80
  )
  gridExtra::grid.arrange(p1, p2, ncol = 1, nrow = 2)
  dev.off()
  rm(p1, p2)
  
  return(dat)
}

normalize_dat_one_T = function(dat_one_Temperature) {
  # store the original unnormalized data in another column
  dat_one_Temperature$raw_log_Result_Value = dat_one_Temperature$log_Result_Value
  # replace log_Result_Value by normalized value
  dat_one_Temperature$log_Result_Value =
    normalizeByMedian(dat_one_Temperature$log_Result_Value,
                      dat_one_Temperature$SampleName)
  return(dat_one_Temperature)
}

normalizeByMedian <- function(logIntensity, sampleGroup) {
  overallMedian <- median(logIntensity, na.rm = TRUE)
  normalized_logIntensity = overallMedian + ave(
    logIntensity,
    sampleGroup,
    FUN = function(x)
      x - median(x, na.rm = TRUE)
  )
  return(normalized_logIntensity)
}

normalize_oneProject_2dCETSA = normalize_oneProject_CETSA
