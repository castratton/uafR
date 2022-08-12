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

dictionate = function(compounds){
  compound_cids = get_cid(compounds)
  compounds_SDF = SDFset()
  if(any(is.na(compound_cids$cid)) == T){
    print("Unmatched Compounds Detected!! Replacing Compounds that have an 'NA' with Acetone.")  #compound_cids, n = length(compound_cids$cid))
    compound_cid_set = toString(compound_cids[[2]])
    compound_cid_set = gsub("NA", "180", compound_cid_set)
    compound_cid_set = gsub("[c\\\"() ]","",compound_cid_set)
    compound_cid_set = gsub("\n", "", compound_cid_set)
    if(length(compound_cids$cid) > 135){
      compound_cid_list = split(compound_cids, ceiling(seq_along(compound_cids[[2]])/135))
      for(i in 1:length(compound_cid_list)){
        compound_cid_set = toString(compound_cid_list[[i]][2])
        compound_cid_set = gsub("NA", "180", compound_cid_set)
        compound_cid_set = gsub("[c\\\"() ]","",compound_cid_set)
        compound_cid_set = gsub("\n", "", compound_cid_set)
        url_compounds = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',compound_cid_set,'/SDF')
        compounds_SDF_tmp = read.SDFset(url_compounds)
        compounds_SDF = append(compounds_SDF, compounds_SDF_tmp)
      }
    } else {
      compound_cid_set = toString(compound_cids[[2]])
      compound_cid_set = gsub("NA", "180", compound_cid_set)
      compound_cid_set = gsub("[c\\\"() ]","",compound_cid_set)
      compound_cid_set = gsub("\n", "", compound_cid_set)
      url_compounds = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',compound_cid_set,'/SDF')
      compounds_SDF = read.SDFset(url_compounds)
    }
    # warning("Unmatched Compounds Detected!! Remove or Replace Compounds that have an 'NA' as their CID!")
  } else {
    compound_cids = compound_cids
    if(length(compound_cids$cid) > 135){
      compound_cid_list = split(compound_cids, ceiling(seq_along(compound_cids[[2]])/135))
      for(i in 1:length(compound_cid_list)){
        compound_cid_set = toString(compound_cid_list[[i]][2])
        compound_cid_set = gsub("NA", "180", compound_cid_set)
        compound_cid_set = gsub("[c\\\"() ]","",compound_cid_set)
        compound_cid_set = gsub("\n", "", compound_cid_set)
        url_compounds = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',compound_cid_set,'/SDF')
        compounds_SDF_tmp = read.SDFset(url_compounds)
        compounds_SDF = append(compounds_SDF, compounds_SDF_tmp)
      }
    } else {
      compound_cid_set = toString(compound_cids[[2]])
      compound_cid_set = gsub("NA", "180", compound_cid_set)
      compound_cid_set = gsub("[c\\\"() ]","",compound_cid_set)
      compound_cid_set = gsub("\n", "", compound_cid_set)
      url_compounds = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',compound_cid_set,'/SDF')
      compounds_SDF = read.SDFset(url_compounds)
    }
  }
  return(compounds_SDF)
}

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

# Looping through every compound and testing against the dictionaries. Will be converted to a function 
# that allows users to define their own dictionaries and Tanimoto thresholds.

TanimotoR = function(SDF_input, comparison_set, target_value){
  output = c()
  for(j in 1:length(SDF_input)){
    batch_test_set = fmcsBatch(SDF_input[[j]],comparison_set)
    # print(batch_test_set)
    if (max(as.vector(batch_test_set[,4])) > target_value)
    {
      output[j] = "Match"
    }else
    {
      output[j] = "No Match"
    }
  }
  return(output)
}

########################################################################################################

standardify = function(u, v, w, x, y, z){
  standardised = (u/v*w)/x/y/z
}

standardifyIt = function(data_in, int_standard, IS_molecules, IS_ng, collect_time, sample_amt){
  data_in_standard = data_in[data_in$cmps_correct != int_standard, c("cmps_correct","mean_RT")]
  for(i in 7:ncol(data_in)){
    Input = data_in[,i][data_in$cmps_correct != int_standard]
    IS = sum(data_in[,i][data_in$cmps_correct == int_standard])
    if(IS > 0){
      tmp = standardify(u = Input, 
                        v = IS, 
                        w = IS_molecules, 
                        x = IS_ng, 
                        y = collect_time, 
                        z = sample_amt)
      data_in_standard = data.frame(data_in_standard, 
                                    tmp)
    } else {
      IS = sum(data_in[,i-1][data_in$cmps_correct == int_standard])
      tmp = standardify(u = Input, 
                        v = IS, 
                        w = IS_molecules, 
                        x = IS_ng, 
                        y = collect_time, 
                        z = sample_amt)
      data_in_standard = data.frame(data_in_standard, 
                                    tmp)
    }
  }
  colnames(data_in_standard) = colnames(data_in[,-c(1:4)])
  return(data_in_standard)
}

#*** PROPRIETARY - DO NOT DISTRIBUTE **