#*** PROPRIETARY - DO NOT DISTRIBUTE **
#Unknowns Analysis Data Processing Script (beta version) 

##Table exported from Unknowns must include 'File Name' or 'Sample Name' info and be saved as a .csv in a known directory
##Required Packages !!

library(data.table)
library(GCalignR)
library(expss)
library(dplyr)
library(tidyr)
library(webchem)
library(ChemmineR)
library(fmcsR)

##Loading the raw Unknowns analysis output
unknowns_all = read.csv("") #This path and filename need to be updated to your input file

out = unknowns_all %>%
  group_split(File.Name) ##Must be changed to 'Sample.Name' if that is the preferred identifier !!

# This is no longer necessary for the alignment to occur, but users may still want to have it as an option (uaf2csv)
# for (i in 1:length(out)){
#   write.csv(out[[i]],paste0("C:/Users/cstratton/Desktop/2019_The_Land_Institute/GCMS_Analyzed/Kernza_Intercropping/CSV_data/Alignment_Split/",unique(out[[i]]$File.Name),".csv"))
# } ## User must manually create a new folder where the individual .csv files for each sample will be saved

fileNames = as.character(unique(unknowns_all$File.Name))
fileNames = gsub('\\.','_', fileNames)
fileNames = gsub('-','_', fileNames)
fileNames = gsub('\\+','plus', fileNames)
#Rename list items of GCalignment input file
names(out) <- fileNames

# CHecking input for alignment
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
cmp_matrix = aligned_peaks$aligned$Compound.Name[,-1]
probs_matrix = aligned_peaks$aligned$Match.Factor[,-1]

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

#Communicating with PubChem to bring in SDF files for every identified/published compound in every sample
## Download/URL restrictions!!!! --> length of URL and time it takes to download data limits
SDF_list = list()

for(i in 1:length(out)){
  cids_tmp = get_cid(out[[i]][[5]])
  cid_set = toString(cids_tmp[[2]])
  cid_set = gsub(" ","",cid_set)
  cid_set = gsub("NA","180",cid_set)
  url_tmp = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',cid_set,'/SDF')
  SDF_set = read.SDFset(url_tmp)
  SDF_list[i] = SDF_set
}
names(SDF_list) = fileNames

# What the compound dictionary inputs look like [dictionate()]
contaminants = c("Perfluorotributylamine","Methanol",
                 "Acetone","Benzene",
                 "Toluene","M-Xylene",
                 "O-Xylene","P-Xylene",
                 "Trichloroethane","dimethoxy(dimethyl)silane")

GLVs = c("(Z)-3-hexenal", "(E)-2-hexenal", #Green Leafy Volatiles
         "(Z)-3-hexen-1-ol", "(Z)-3-hexen-1-yl acetate")

acAs = c("Vinblastine",  #Anti-Cancer Alkaloids - "Vincrisine",
         "Vinorelbine","Vindesine",
         "Taxol","Rohitukine",
         "Homoharringtonine","Ellipticine",
         "Acronycine","Camptothecin",
         "Thalicarpine")

Alkaloids = c("Trimethylamine","Dimethylamine")

# Acquiring SDFs for compound dictionaries [dictionate()]
## Contaminants
contaminant_cids = get_cid(contaminants)
contaminant_cid_set = toString(contaminant_cids[[2]])
contaminant_cid_set = gsub(" ","",contaminant_cid_set)
url_contaminants = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',contaminant_cid_set,'/SDF')
contaminant_SDF_set = read.SDFset(url_contaminants)

## Green Leafy Volatiles
GLV_cids = get_cid(GLVs)
GLV_cid_set = toString(GLV_cids[[2]])
GLV_cid_set = gsub(" ","",GLV_cid_set)
url_GLVs = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',GLV_cid_set,'/SDF')
GLV_SDF_set = read.SDFset(url_GLVs)

## Anti-Cancer alkaloids
acA_cids = get_cid(acAs)
acA_cid_set = toString(acA_cids[[2]])
acA_cid_set = gsub(" ","",acA_cid_set)
url_acAs = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',acA_cid_set,'/SDF')
acA_SDF_set = read.SDFset(url_acAs)

# Alkaloids
Alkaloid_cids = get_cid(Alkaloids)
Alkaloid_cid_set = toString(Alkaloid_cids[[2]])
Alkaloid_cid_set = gsub(" ","",Alkaloid_cid_set)
url_Alkaloids = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',Alkaloid_cid_set,'/SDF')
Alkaloid_SDF_set = read.SDFset(url_Alkaloids)

# Adding new columns to the split datasets to fill with results from dictionary matching.
out2 = Map(function(x,y){
  x$Contaminant = y 
  x$GLV = y 
  x$acA = y 
  x$Alkaloid = y 
  x}, 
  out, 
  "Empty")

# Looping through every compound and testing against the dictionaries. Will be converted to a function 
# that allows users to define their own dictionaries and Tanimoto thresholds.
for(i in 1:length(SDF_list)){
  current_set = SDF_list[[i]] 
  for(j in 1:length(current_set)){
    batch_test_set = fmcsBatch(current_set[[j]],contaminant_SDF_set) 
    batch_test_set_GLV = fmcsBatch(current_set[[j]],GLV_SDF_set) 
    batch_test_set_acA = fmcsBatch(current_set[[j]],acA_SDF_set) 
    batch_test_set_Alkaloid = fmcsBatch(current_set[[j]],Alkaloid_SDF_set) 
    if (max(as.vector(batch_test_set[,4])) > 0.999){
      out2[[i]][j,ncol(out[[i]])+1] = "Contaminant"} 
    else {
      out2[[i]][j,ncol(out[[i]])+1] = "Okay"} 
    if (max(as.vector(batch_test_set_GLV[,4])) > 0.85){
      out2[[i]][j,10] = "GLV"} 
    else {
      out2[[i]][j,10] = "No"} 
    if (max(as.vector(batch_test_set_acA[,4])) > 0.75){
      out2[[i]][j,11] = "acA"} 
    else {out2[[i]][j,11] = "No"} 
    if (max(as.vector(batch_test_set_Alkaloid[,4])) > 0.97){
      out2[[i]][j,12] = "Alkaloid"} 
    else {out2[[i]][j,12] = "No"}
  }
  }

#*** PROPRIETARY - DO NOT DISTRIBUTE **