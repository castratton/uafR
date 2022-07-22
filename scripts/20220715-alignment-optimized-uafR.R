#*** PROPRIETARY - DO NOT DISTRIBUTE **
#Unknowns Analysis Data Processing Script (beta version) 

##Table exported from Unknowns must include 'File Name' or 'Sample Name' info and be saved as a .csv in a known directory
##Required Packages !!

library(data.table)
library(GCalignR)
library(expss)
library(dplyr)
library(tidyr)

##Loading the raw Unknowns analysis output
unknowns_all = read.csv("C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Morrison_Demo/CSV_data/Aggregated_All/") #This path and filename need to be updated to your input file

out = unknowns_all %>%
  group_split(File.Name) ##Must be changed to 'Sample.Name' if that is the preferred identifier !!

# This is no longer necessary for the alignment to occur, but users may still want to have it as an option
# for (i in 1:length(out)){
#   write.csv(out[[i]],paste0("C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Kernza_Intercropping/CSV_data/Alignment_Split/",unique(out[[i]]$File.Name),".csv"))
# } ## User must manually create a new folder where the individual .csv files for each sample will be saved

# Importing dataset into R. Creating a list of file names, then reading the relevant columns into GC_input
# fileList = list.files(path="C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Kernza_Intercropping/CSV_data/Alignment_Split/", pattern=".csv", full.names = TRUE)
# GC_input = lapply(fileList, fread, select = c("Component.RT", "Component.Area"))

# GC_input is without list item labels, the next step imports labels from the file names
# This code can probably be optimized

# fileNames = list.files(path="C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Kernza_Intercropping/CSV_data/Alignment_Split/", 
#                        pattern = NULL, all.files = TRUE, 
#                        full.names = FALSE, recursive = FALSE,
#                        ignore.case = FALSE, include.dirs = FALSE, 
#                        no.. = FALSE)
# 
# fileNames = substr(fileNames,1,nchar(fileNames)-4)
# fileNames = data.frame(fileNames)
# fileNames = fileNames[-c(1:2),]
# fileNames = gsub('\\.','_', fileNames)
# fileNames = gsub('-','_', fileNames)

fileNames = as.character(unique(unknowns_all$File.Name))

# fileNames = substr(fileNames,1,nchar(fileNames)-4)
# fileNames = data.frame(fileNames)
# fileNames = fileNames[-c(1:2),]
fileNames = gsub('\\.','_', fileNames)
fileNames = gsub('-','_', fileNames)

#Rename list items of GCalignment input file
# names(GC_input) <- fileNames
names(out) <- fileNames

#Rename list items of GCalignment input file
# names(GC_input) <- fileNames
# names(out) <- fileNames
# GC_input=GC_input[-127]
check_input(out)


##Run GCMS alignment

#Determine appropriate value for min_diff_peak2peak
tmp = peak_interspace(data = out, 
                      rt_col_name = "Component.RT", 
                      quantile_range = c(0, 0.95),
                      by_sample = F)  #Could potentially extract a value for each sample individually
shift = tmp$Summary[4]                     #{
peak2mean = tmp$Summary[6]-tmp$Summary[4]  #{ <- Current setup is working well enough on our data
peak2peak_min = tmp$Summary[2]             #{

#Run alignment <-- LOTS OF ROOM FOR OPTIMIZATION (data science "experiment" - 
# test a regular series of values here to find best fit values across data sets)
aligned_peaks<- align_chromatograms(out, 
                                    rt_col_name = "Component.RT", 
                                    max_linear_shift = shift, 
                                    # max_diff_peak2mean = peak2mean, 
                                    min_diff_peak2peak = peak2peak_min)

#Check output
gc_heatmap(aligned_peaks)
plot(aligned_peaks)
print(aligned_peaks)
View(aligned_peaks$aligned)
#Write output to .csv
# write.csv(aligned_peaks$aligned, 
#           "C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Alfalfa_Intercropping/20220531-alfalfa-aligned.csv")

options(scipen = 999)
#Adding compound names and probabilities
# Cmpd_lookup = lapply(fileList, fread) %>% rbindlist()
#Cmpd_lookup <- Cmpd_lookup %>%
# unite("Compound; Prob", Compound:Probability, sep = "; ")
# Cmpd_lookup = Cmpd_lookup[,c(2,3,4)]
# cmp_tmp = as.data.frame(Cmpd_lookup[,c("Component.RT","Component.Area","Compound.Name","Match.Factor")])
cmp_tmp = as.data.frame(unknowns_all[,c("Component.RT","Component.Area","Compound.Name","Match.Factor")])
#For Excel method: .csv file
# write.csv(Cmpd_lookup, "C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Alfalfa_Intercropping/GCAligned/20220531-alfalfa-CMP-list.csv")
#write.csv(Cmpd_lookup, "C:/Intercropping/Output/RTCmpds.csv")

Orig=aligned_peaks$aligned$Component.Area
# vlookup(Orig[1,3],cmp_tmp, result_column = 2, lookup_column = 1)
orig_new = Orig[,-1]

ret_time = nrow(orig_new)
sample = ncol(orig_new)

cmp_matrix = data.frame(matrix(ncol = sample, 
                               nrow = ret_time))
probs_matrix = data.frame(matrix(ncol = sample, 
                                 nrow = ret_time))

for (i in 1:ret_time){
  for (j in 1:sample){
    cmp_matrix[i,j] = as.character(vlookup(orig_new[i,j], 
                                           cmp_tmp, 
                                           result_column = 3, 
                                           lookup_column = 2))
    probs_matrix[i,j] = vlookup(orig_new[i,j], 
                                cmp_tmp, 
                                result_column = 4, 
                                lookup_column = 2)
  }
}

# cbind(t(cmp_matrix[1,]),t(probs_matrix[1,]))
best_cmps = vector()

for (i in 1:nrow(probs_matrix)){
  rt_set = as.data.frame(cbind(t(probs_matrix[i,]),t(cmp_matrix[i,])))
  row.names(rt_set) = NULL
  colnames(rt_set) = c("Probs","Compounds")
  ordered_set = rt_set[order(rt_set$Probs),]
  ordered_set = ordered_set[complete.cases(ordered_set),]
  max_tmp = tail(ordered_set, n = 1)
  best_cmps[i] = as.character(max_tmp[[2]])
}
cmps_correct = unlist(best_cmps)

aligned_with_CMPS = cbind(cmps_correct, Orig)
View(aligned_with_CMPS)

write.csv(aligned_with_CMPS,"C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Silphium_InterchemHerbivore/20220718-interchemHerbivore-ALL-R-CMPs.csv")
##########################################################################################################
# ************* SHOULD SAVE FILE HERE !!!!!!!!!!!!!!! ************ #######################################
##########################################################################################################

#Additional steps we could take:
##-Merge duplicates <-- need to average retention times and sum areas
##-Remove contaminants <-- m/z value matching for most common (septa bleed etc.)
##-Quantify relative quantities (relative to internal standard <-- Rob feedback required)
##-Summarize based on key sub-structures (e.g alkaloids, terpenoids, aromatics, etc.) <-- could develop "dictionaries" for specific fields
##--Would be an optional argument in a function. Users would have options to select (e.g. "plant", "medical", "others...")

#*** PROPRIETARY - DO NOT DISTRIBUTE **