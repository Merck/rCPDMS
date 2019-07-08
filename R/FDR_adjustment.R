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



## FDR adjustment functions

# Pulldown ----------------------------------------------------------------

add_FDR_pulldown = function(results_input, baseline = 0){
  if(!is.data.frame(results_input)){
    cat("================= Adding FDR correction to pVal for non-baseline results ================= \n")
    results_list = lapply(results_input, add_FDR_pulldown, baseline)
    return(results_list)
  }else{
    this_project = unique(as.character(results_input$Project_Experiment_SampleGroup))
    cat("+++ Project_Experiment_SampleGroup: ",this_project," +++ \n")
    results_all = add_FDR_one_SampleGroup_pulldown(results_input)
    return(results_all)
  }
}

add_FDR_one_SampleGroup_pulldown = function(results_SampleGroup){
  pval_cols = names(results_SampleGroup)[grepl("pVal_", names(results_SampleGroup))]
  pval_cols_fdr = sapply(pval_cols, function(x) paste0(x,"_FDR"))
  for(i in 1:length(pval_cols)){
    results_SampleGroup[[pval_cols_fdr[i]]] = p.adjust(results_SampleGroup[[pval_cols[i]]], method = "BH")
  }
  return(results_SampleGroup)
}



# Pro10K ------------------------------------------------------------------

add_FDR_pro10K = function(results_input, baseline = "DMSO"){
  if(!is.data.frame(results_input)){
    cat("================= Adding FDR correction to pVal for non-baseline results ================= \n")
    results_list = lapply(results_input, add_FDR_pro10K, baseline)
    return(results_list)
  }else{
    this_project = unique(as.character(results_input$Project_Experiment))
    cat("+++ Project_Experiment: ",this_project," +++ \n")
    
    results_split = split(results_input, results_input$SampleGroup)
    results_all_list = lapply(results_split, add_FDR_one_SampleGroup_pro10K, baseline)
    results_all = bind_rows(results_all_list)
    return(results_all)
  }
}

add_FDR_one_SampleGroup_pro10K = function(results_SampleGroup, baseline = "DMSO"){
  if (unique(results_SampleGroup$SampleGroup) == baseline) {
    results_SampleGroup$pVal_FDR <- results_SampleGroup$pVal_limma_FDR <- NA
  } else {
    results_SampleGroup$pVal_FDR = p.adjust(results_SampleGroup$pVal, method = "BH")
    results_SampleGroup$pVal_limma_FDR = p.adjust(results_SampleGroup$pVal_limma, method = "BH")
  }
  return(results_SampleGroup)
}

# add_FDR_one_SampleGroup_pro10K = function(results_SampleGroup, baseline = "DMSO"){
#   df = results_SampleGroup %>% select(Project_Experiment,proteinID, SampleGroup, pVal) %>% 
#     filter(SampleGroup !=baseline & !is.na(pVal)) %>% 
#     arrange(pVal) 
#   if(nrow(df)==0){
#     results_SampleGroup$pVal_FDR = NA
#     return(results_SampleGroup)
#   } 
#   df$pVal_order = 1:nrow(df)
#   df$pVal_FDR = df$pVal * nrow(df) / df$pVal_order
#   df = df %>% select(-pVal, -pVal_order)
#   
#   results_SampleGroup = results_SampleGroup %>%
#     left_join(df, by = c("Project_Experiment","proteinID", "SampleGroup"))
#   return(results_SampleGroup)
# }


# CETSA -------------------------------------------------------------------

add_FDR_CETSA = function(results_input, baseline = "DMSO"){
  if(!is.data.frame(results_input)){
    cat("================= Adding FDR correction to pVal for non-baseline results ================= \n")
    results_list = lapply(results_input, add_FDR_CETSA, baseline)
    return(results_list)
  }else{
    this_project = unique(as.character(results_input$Project_Experiment))
    cat("+++ Project_Experiment: ",this_project," +++ \n")
    
    results_split = split(results_input, results_input$SampleGroup)
    results_all_list = lapply(results_split, add_FDR_one_SampleGroup_CETSA, baseline)
    results_all = bind_rows(results_all_list)
    return(results_all)
  }
}

add_FDR_one_SampleGroup_CETSA = function(results_SampleGroup, baseline = "DMSO") {
  if(unique(results_SampleGroup$SampleGroup) == baseline){
    results_SampleGroup$pVal_IC50_diff_FDR = NA
  }else{
    results_SampleGroup$pVal_IC50_diff_FDR = p.adjust(results_SampleGroup$pValue_IC50_diff, method = "BH")
  }
  return(results_SampleGroup)
}

add_FDR_2dCETSA = add_FDR_CETSA
add_FDR_2dCETSA_by_temp = add_FDR_pulldown
