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



rm(list = ls())
library(rCPDMS)

# pulldown ----------------------------------------------------------------
workflow = "pulldown" 
dataPath = "../data_raw/pulldown_Sample.txt"
saveFolder = file.path("output",workflow)
baseline = 0

outlier = FALSE
resultType = "Intensity" 
normalization = FALSE
imputation = FALSE
results = MSanalysis(workflow, dataPath, resultType, baseline, saveFolder, outlier, normalization, imputation)

# CETSA -------------------------------------------------------------------
workflow = "CETSA" 
dataPath = "../data_raw/CETSA_Sample.txt"
saveFolder = file.path("output",workflow)
baseline = "DMSO" 

outlier = FALSE
resultType = "Intensity" 
normalization = TRUE
imputation = TRUE
results = MSanalysis(workflow, dataPath, resultType, baseline, saveFolder, outlier, normalization, imputation)

# 2d CETSA ----------------------------------------------------------------
workflow = "2dCETSA" 
dataPath = "../data_raw/2dCETSA_Sample.txt"
saveFolder = file.path("output",workflow)
baseline = "DMSO" 

outlier = FALSE
resultType = "Intensity" 
normalization = FALSE
imputation = TRUE
results = MSanalysis(workflow, dataPath, resultType, baseline, saveFolder, outlier, normalization, imputation)
