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



#########################################################################
#
# The code will take the maxquant protein output and mapping file to generate
# input data file for CPDMS package. 
#
#
##########################################################################


require(reshape)

proteinfile <- "protein.txt"
mappingfile <- "MappingTemplaete.csv"
outfile <- "test_rawdata.txt"

map.dat <- read.csv(mappingfile, header=T)
#Use this command to remove columns that are entirely NA values, it will elave columns where only some vlaues are NA
map.dat <- map.dat[ , ! apply(map.dat, 2 , function(x) all(is.na(x)) ) ]
map.dat["Result_Column_Name"] <- t(sub(" ", ".", t(map.dat["Result_Column_Name"] )))
protein.dat <- read.table(proteinfile, header=T, sep='\t')
protein.dat <- data.frame(protein.dat)


col.keep <- c("Majority.protein.IDs", sub(" ", ".", t(map.dat["Result_Column_Name"])))
protein.dat <- protein.dat[which(protein.dat["Majority.protein.IDs"] != ""), col.keep]
names(protein.dat)[1] <-"ProteinID"

protein.dat.unpivot <- melt(protein.dat, id = c("ProteinID"))
names(protein.dat.unpivot) <- c("ProteinID", "Result_Column_Name", "Result_Value")

protein.data.new <- merge(protein.dat.unpivot, map.dat, by.x="Result_Column_Name", by.y="Result_Column_Name", all.x=T, all.y=F)
protein.data.new <- data.frame(protein.data.new)

write.table(protein.data.new, outfile, col.names = TRUE, row.names = FALSE, sep = '\t')
