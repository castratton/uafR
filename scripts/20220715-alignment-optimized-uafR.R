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
library(ggplot2)
library(jsonlite)

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

## ***************************************************************************************************** ##
##     mz_splitter() - separates overlapping peaks and accurately assigns identity to partial matches    ##
##     ***Could use some work still but basically there. Still need to incorporate TanimotoMerger***     ##
## ***************************************************************************************************** ##
# input = aligned_peaks
# samples_with = 3

mz_splitter = function(input, samples_with){
  floor_dec = function(x, level=1) round(x - 5*10^(-level-1), level)
  
  Orig=input$aligned$Component.Area

  orig_new = Orig[,-1]
  
  ret_time = nrow(orig_new)
  sample_name = ncol(orig_new)
  
  new_row = rep(0, sample_name)
  rt_row = mz_row = probs_row = cmp_row = new_row
  rt_row2 = mz_row2 = probs_row2 = cmp_row2 = new_row
  
  cmp_matrix = as.data.frame(matrix(unlist(input$aligned$Compound.Name[,-1]), 
                                    nrow = ret_time, 
                                    byrow = F,
                                    dimnames = list(input$aligned$Compound.Name[,1], 
                                                    colnames(input$aligned$Compound.Name[,-1]))), 
                             stringsAsFactors = F)
  for (cmp_col in 1:ncol(cmp_matrix)){
    if (is.na(cmp_matrix[,cmp_col][is.na(suppressWarnings(as.integer(cmp_matrix[,cmp_col])))])[5]){
      cmp_matrix[,cmp_col] = input$aligned$Compound.Name[,cmp_col+1]
    } else {
    }
    cmp_matrix = as.data.frame(cmp_matrix)
  }
  
  
  probs_matrix = as.data.frame(matrix(unlist(input$aligned$Match.Factor[,-1]), 
                                      nrow = ret_time, 
                                      byrow = F, 
                                      dimnames = list(input$aligned$Compound.Name[,1],
                                                      colnames(input$aligned$Compound.Name[,-1]))), 
                               stringsAsFactors = F)
  
  mz_matrix = as.data.frame(matrix(unlist(input$aligned$Base.Peak.MZ[,-1]), 
                                   nrow = ret_time, 
                                   byrow = F, 
                                   dimnames = list(input$aligned$Compound.Name[,1],
                                                   colnames(input$aligned$Compound.Name[,-1]))), 
                            stringsAsFactors = F)
  
  rt_matrix = as.data.frame(matrix(unlist(input$aligned$Component.RT[,-1]), 
                                   nrow = ret_time, 
                                   byrow = F, 
                                   dimnames = list(input$aligned$Compound.Name[,1],
                                                   colnames(input$aligned$Compound.Name[,-1]))), 
                            stringsAsFactors = F)
  
  area_matrix = as.data.frame(matrix(unlist(input$aligned$Component.Area[,-1]), 
                                     nrow = ret_time, 
                                     byrow = F, 
                                     dimnames = list(input$aligned$Compound.Name[,1],
                                                     colnames(input$aligned$Compound.Name[,-1]))), 
                              stringsAsFactors = F)
  
  cmp_matrix_multiples = probs_matrix_multiples = mz_matrix_multiples = rt_matrix_multiples = area_matrix_multiples = NULL
  cmp_matrix_solos = probs_matrix_solos = mz_matrix_solos = rt_matrix_solos = area_matrix_solos = NULL
  cmp_matrix_none = probs_matrix_none = mz_matrix_none = rt_matrix_none = area_matrix_none = NULL

  rt_column = c()

  for (i in 1:ret_time){
    mz_current = unique(unlist(mz_matrix[i,]))
    probs_current = unique(unlist(probs_matrix[i,]))
    cmp_current = unique(cmp_matrix[i, , drop = T])
    cmp_names = rep("NA", length(cmp_current))
    for (cmp in 1:length(cmp_current)){
      cmp_tmp = cmp_current[[cmp]]
      cmp_names[cmp] = as.character(cmp_tmp)
    }
    cmp_names = cmp_names[cmp_names!="0"]
    cmp_current = cmp_names[!is.na(cmp_names)]
    rt_current = unique(unlist(rt_matrix[i,]))
    rt_current = unique(floor_dec(rt_current, 2))
    area_current = unique(unlist(area_matrix[i,]))
    area_current = unique(floor_dec(area_current, 2))
    
    count_mz = c()
    for (k in 1:length(mz_current)){
      count_tmp = sum(as.numeric(mz_matrix[i,]==mz_current[k]))
      count_mz = c(count_mz, 
                   count_tmp)
    }
    mz_counts = as.data.frame(cbind(count_mz, 
                                    as.numeric(unlist(mz_current))))
    colnames(mz_counts) = c("count", 
                            "mz_current")
    
    count_probs = c()
    for (k in 1:length(probs_current)){
      count_tmp = sum(as.numeric(probs_matrix[i,]==probs_current[k]))
      count_probs = c(count_probs, 
                      count_tmp)
    }
    probs_counts = as.data.frame(cbind(count_probs, 
                                       as.numeric(unlist(probs_current))))
    colnames(probs_counts) = c("count", 
                               "probs_current")
    
    count_cmp = c()
    for (k in 1:length(cmp_current)){
      cmp_going = cmp_matrix[i, , drop = T]
      cmp_names = rep("NA", length(cmp_going))
      for (cmp in 1:length(cmp_going)){
        cmp_tmp = cmp_going[[cmp]]
        cmp_names[cmp] = as.character(cmp_tmp)
      }
      cmp_names = cmp_names[cmp_names!="0"]
      cmp_going = cmp_names[!is.na(cmp_names)]
      count_tmp = sum(as.numeric(cmp_going==cmp_current[k]))
      count_cmp = c(count_cmp, 
                    count_tmp)
    }
    cmp_counts = as.data.frame(cbind(count_cmp, 
                                     cmp_current))
    colnames(cmp_counts) = c("count", 
                             "cmp_current")
    
    count_rt = c()
    for (k in 1:length(rt_current)){
      count_tmp = sum(floor_dec(as.numeric(as.vector(rt_matrix[i,])), 
                                2)==rt_current[k])
      count_rt = c(count_rt, 
                   count_tmp)
    }
    rt_counts = as.data.frame(cbind(count_rt,as.numeric(unlist(rt_current))))
    colnames(rt_counts) = c("count", 
                            "rt_current")
    
    count_area = c()
    for (k in 1:length(area_current)){
      count_tmp = sum(floor_dec(as.numeric(as.vector(area_matrix[i,])), 
                                2)==as.numeric(area_current[k]))
      count_area = c(count_area, 
                     count_tmp)
    }
    area_counts = as.data.frame(cbind(count_area,as.numeric(unlist(area_current))))
    colnames(area_counts) = c("count", 
                              "area_current")
    
    mz_multiples = mz_counts$mz_current[which(mz_counts$count >= samples_with & mz_counts$mz_current != 0)]
    
    cmp_multiNew = cmp_matrix[i, , drop = T]
    cmp_multiNames = rep("NA", length(cmp_multiNew))
    for (cmp in 1:length(cmp_multiNew)){
      cmp_tmp = cmp_multiNew[[cmp]]
      cmp_multiNames[cmp] = as.character(cmp_tmp)
    }
    multiples_identities = as.data.frame(cbind(cmp_multiNames, 
                                               unlist(probs_matrix[i,]),
                                               unlist(mz_matrix[i,]),
                                               unlist(rt_matrix[i,]),
                                               unlist(area_matrix[i,])))
    colnames(multiples_identities) = c("Compound",
                                       "Prob",
                                       "mz",
                                       "RT",
                                       "Area")
    probs_multiples = multiples_identities$Prob[multiples_identities$mz %in% mz_multiples]
    cmp_multiples = multiples_identities$Compound[multiples_identities$mz %in% mz_multiples]
    rt_multiples = multiples_identities$RT[multiples_identities$mz %in% mz_multiples]
    area_multiples = multiples_identities$Area[multiples_identities$mz %in% mz_multiples]
    
    multiples_logical = is.null(mz_multiples)
    multiples_logical2 = !is.na(mz_multiples[1])
    
    mz_solos = mz_counts$mz_current[which(mz_counts$count < samples_with & mz_counts$mz_current != 0)]
    
    solos_identities = as.data.frame(cbind(cmp_multiNames, 
                                           unlist(probs_matrix[i,]),
                                           unlist(mz_matrix[i,]),
                                           unlist(rt_matrix[i,]),
                                           unlist(area_matrix[i,])))
    colnames(solos_identities) = c("Compound",
                                   "Prob",
                                   "mz",
                                   "RT",
                                   "Area")
    probs_solos = solos_identities$Prob[solos_identities$mz %in% mz_solos]
    cmp_solos = solos_identities$Compound[solos_identities$mz %in% mz_solos]
    rt_solos = solos_identities$RT[solos_identities$mz %in% mz_solos]
    area_solos = solos_identities$Area[solos_identities$mz %in% mz_solos]

    solos_logical = is.null(mz_solos)
    solos_logical2 = !is.na(mz_solos[1])
    if (!multiples_logical & multiples_logical2){
      for (l in 1:length(mz_multiples)){
        multiples_matrix = as.data.frame(cbind(cmp_multiNames, 
                                               unlist(probs_matrix[i,]),
                                               unlist(mz_matrix[i,]),
                                               unlist(rt_matrix[i,]),
                                               unlist(area_matrix[i,])))
        colnames(multiples_matrix) = c("Compound",
                                       "Prob",
                                       "mz",
                                       "RT",
                                       "Area")
        
        mz_row = unlist(mz_matrix[i,])
        probs_row = unlist(probs_matrix[i,])
        cmp_row = cmp_multiNames
        rt_row = unlist(rt_matrix[i,])
        area_row = unlist(area_matrix[i,])
        
        mz_rows = probs_rows = cmp_rows = rt_rows = area_rows = rep("NA", length(mz_row))
        # l = 1
        area_yes = rt_yes = cmp_yes = probs_yes = mz_yes = multiples_matrix$mz == mz_multiples[l]
        area_no = rt_no = cmp_no = probs_no = mz_no = multiples_matrix$mz != mz_multiples[l]
        
        mz_row[mz_no] <- 0
        mz_row[mz_yes] <- mz_multiples[l]
        
        max_location = multiples_matrix$mz == mz_multiples[l]
        max_prob = max(as.numeric(as.vector(unlist(multiples_matrix$Prob[max_location]))))
        probs_row[probs_no] <- 0
        probs_row[probs_yes] <- max_prob  
        
        best_cmp = multiples_matrix$Compound[multiples_matrix$Prob == max_prob]
        cmp_row[cmp_no] <- "NA"
        cmp_row[cmp_yes] <- as.character(unlist(best_cmp))
        
        best_rt = multiples_matrix$RT[multiples_matrix$Prob == max_prob]
        rt_row[rt_no] <- 0
        rt_row[rt_yes] <- as.numeric(as.vector(unlist(best_rt)))
        
        # best_area = multiples_matrix$Area[multiples_matrix$Prob == max_prob] ### *** Area needs to be fixed!! Unique values have to be assigned.
        area_row[area_no] <- 0
        # area_row[area_yes] <- as.numeric(as.vector(unlist(best_area)))
        cmp_rows = rbind(cmp_rows, cmp_row)
        probs_rows = rbind(probs_rows, probs_row)
        mz_rows = rbind(mz_rows, mz_row)
        rt_rows = rbind(rt_rows, rt_row)
        area_rows = rbind(area_rows, area_row)
      }
      cmp_matrix_multiples = rbind(cmp_matrix_multiples, 
                                   cmp_rows)
      probs_matrix_multiples = rbind(probs_matrix_multiples, 
                                     probs_rows)
      mz_matrix_multiples = rbind(mz_matrix_multiples, 
                                  mz_rows)
      rt_matrix_multiples = rbind(rt_matrix_multiples, 
                                  rt_rows)
      area_matrix_multiples = rbind(area_matrix_multiples, 
                                    area_rows)} else 
                                      if (!solos_logical & solos_logical2){
                                        solos_matrix = as.data.frame(cbind(cmp_multiNames, 
                                                                           unlist(probs_matrix[i,]),
                                                                           unlist(mz_matrix[i,]),
                                                                           unlist(rt_matrix[i,]),
                                                                           unlist(area_matrix[i,])))
                                        colnames(solos_matrix) = c("Compound",
                                                                   "Prob",
                                                                   "mz",
                                                                   "RT",
                                                                   "Area")
                                        
                                        mz_row = unlist(mz_matrix[i,])
                                        probs_row = unlist(probs_matrix[i,])
                                        cmp_row = cmp_multiNames
                                        rt_row = unlist(rt_matrix[i,])
                                        area_row = unlist(area_matrix[i,])
                                        
                                        rt_yes = cmp_yes = probs_yes = mz_yes = solos_matrix$mz %in% mz_solos
                                        area_no = rt_no = cmp_no = probs_no = mz_no = !(solos_matrix$mz %in% mz_solos)
                                        
                                        mz_row[mz_no] <- 0
                                        mz_row[mz_yes] <- mz_multiples[l]
                                        
                                        max_location = solos_matrix$mz %in% mz_solos
                                        max_prob = max(as.numeric(as.vector(unlist(solos_matrix$Prob[max_location]))))
                                        probs_row[probs_no] <- 0
                                        probs_row[probs_yes] <- max_prob  
                                        
                                        best_cmp = solos_matrix$Compound[solos_matrix$Prob == max_prob]
                                        cmp_row[cmp_no] <- "NA"
                                        cmp_row[cmp_yes] <- as.character(unlist(best_cmp))
                                        
                                        best_rt = solos_matrix$RT[solos_matrix$Prob == max_prob]
                                        rt_row[rt_no] <- 0
                                        rt_row[rt_yes] <- as.numeric(as.vector(unlist(best_rt)))
                                        
                                        # best_area = solos_matrix$Area[solos_matrix$Prob == max_prob]
                                        area_row[area_no] <- 0
                                        # area_row[area_yes] <- as.numeric(as.vector(unlist(best_area)))
                                        
                                        cmp_matrix_solos = rbind(cmp_matrix_solos, 
                                                                 cmp_row)
                                        probs_matrix_solos = rbind(probs_matrix_solos, 
                                                                   probs_row)
                                        mz_matrix_solos = rbind(mz_matrix_solos, 
                                                                mz_row)
                                        rt_matrix_solos = rbind(rt_matrix_solos, 
                                                                rt_row)
                                        area_matrix_solos = rbind(area_matrix_solos, 
                                                                  area_row)
                                      } else {
                                        
                                      }
  }
  mz_splitted = list(rbind(cmp_matrix_multiples, cmp_matrix_solos), 
                     rbind(probs_matrix_multiples, probs_matrix_solos), 
                     rbind(mz_matrix_multiples, mz_matrix_solos), 
                     rbind(rt_matrix_multiples, rt_matrix_solos),
                     rbind(area_matrix_multiples, area_matrix_solos))
  names(mz_splitted) = c("Compounds", "Probabilities", "mz value", "RT", "Component Area")
  return(mz_splitted)
}
################################################################################################

########################################################################

## Ready to become a function or merge with existing
probs_matrix = aligned_peaks$aligned$Match.Factor
cmp_matrix = aligned_peaks$aligned$Compound.Name

probs_matrix = mz_splits$Probabilities
cmp_matrix = mz_splits$Compounds
RT_matrix = mz_splits$RT

best_cmps = vector()
mean_RT = vector()
# i = 24
for (i in 1:nrow(probs_matrix)){
  rt_set = as.data.frame(cbind(as.vector(unlist(probs_matrix[i,])),as.vector(unlist(cmp_matrix[i,])),as.vector(unlist(RT_matrix[i,]))))
  row.names(rt_set) = NULL
  colnames(rt_set) = c("Probs","Compounds","RT")
  ordered_set = rt_set[order(rt_set$Probs),]
  ordered_set = ordered_set[complete.cases(ordered_set),]
  max_tmp = tail(ordered_set, n = 1)
  best_cmps[i] = as.character(max_tmp[[2]])
  mean_RT[i] = droplevels(max_tmp[3])
}
cmps_correct = unlist(best_cmps)

aligned_with_CMPS = cbind(cmps_correct, Orig)
##########################################################################################################
# ************* SHOULD SAVE FILE HERE !!!!!!!!!!!!!!! ************ #######################################
##########################################################################################################
contaminants = c("Perfluorotributylamine","Methanol","Benzene",
                 "Toluene","M-Xylene",
                 "O-Xylene","P-Xylene",
                 "Trichloroethane","dimethoxy(dimethyl)silane")

insecticides = c("abamectin", "acephate","acetamiprid",
                 "allethrin","arsenic trioxide","azadirachtin",
                 "bifenthrin","borate",
                 "carbaryl","chlorantraniliprole","clothianidin",
                 "cryolite","cyfluthrin","cypermethrin",
                 "diflubenzuron","dinotefuran","disulfoton",
                 "emamectin B1a","fipronil","fluvalinate",
                 "hydramethylnon","imidacloprid",
                 "lambda-cyhalothrin","malathion","salannin",
                 "permethrin","pyrethrin","resmethrin",
                 "rotenone","silica","spinosad",
                 "sulfluramid","thiamethoxam")

herbicides = c("2,4-D","benefin","bensulide",
               "bentazon","bromoxynil","cacodylic acid",
               "calcium acid methanearsonate","carfentrazone","chlorsulfuron",
               "clethodim","CHEMBL2372957","dicamba",
               "dichlobenil","dimethenamid-P",
               "diquat","dithiopyr","EPTC",
               "fluazifop","fluroxypyr","foramsulfuron",
               "glufosinate","glyphosate","halosulfuron",
               "hexazinone","imazapyr","isoxaben",
               "MCPA","mecoprop","metolachlor",
               "MSMA","napropamide","oryzalin",
               "oxadiazon","oxyfluorfen","pelargonic acid",
               "pendimethalin","prodiamine","pronamide",
               "quinclorac","sethoxydim","siduron",
               "sulfosulfuron","tebuthiuron","triclopyr",
               "trifloxysulfuron-sodium","trifluralin")

pheromones = c("(E)-11-hexadecenal", "(E,E)-10,12-hexadecadienal","(E)-10-hexadecenal", 
               "(Z)-10-hexadecenal", "(E,E)-10,12-hexadecadienal", 
               "(Z)-11-hexadecenyl acetate", "(Z)-11-hexadecenal", 
               "(Z)-11-hexadecenol", "(4aS,7S,7aR)-nepetalactone", "(1R,4aS,7S,7aR)-nepetalactol",
               "(E,Z)-3,13-octadecadienyl acetate", "(Z,Z)-3,13-octadecadienyl acetate", 
               "(E)-11-tetradecenyl acetate", "(E,E)-9,11-tetradecadienyl acetate", "(E)-11-tetradecenol",
               "(E)-11-hexadecenyl acetate", "(Z,E)-9,12-tetradecadienyl acetate", 
               "(Z)-9-tetradecenyl acetate", "(Z)-11-hexadecenyl acetate",
               "(Z,E)-9,12-tetradecadienol", "(Z)-9-tetradecenol", "(Z)-11-hexadecenol",
               "(Z)-9-tetradecenyl acetate", "(Z)-9-tetradecenol", "tetradecyl acetate",
               "(R)-delta-heptalactone", "(Z,E)-5,7-dodecadienol", 
               "(Z,E)-5,7-dodecadienyl propionate")

morrison_matches = c("Pivalaldehyde, semicarbazone", "2-Butenal", 
                     "2,4-Dimethyl-1-heptene","4-hydroxy-4-methyl-5-phenoxy-2-Pentanone", 
                     "1-Octen-3-ol", "Butanal", 
                     "4,5-Dichloro-1,3-dioxolan-2-one", "3-Octanone", "Decane", 
                     "Mesitylene", "Benzene, 1,2,4-trimethyl-", "D-Limonene", 
                     "2-Phenylbutan-1-Ol", "1,4-diethyl-2-Methylbenzene", 
                     "1,2-diethyl-3-Methylbenzene", "1,3,8-p-Menthatriene", 
                     "4-[Dichloromethyl]-2-[[2-[1-Methyl-2-Pyrrolidinyl]Ethyl]Amino-6-Trichloromethylpyrimidine", 
                     "Benzene, (2-methyl-1-propenyl)-", "1-Phenyl-1-butene", 
                     "Linalool", "Nonanal", "2-Thiophenecarboxylic acid, 5-nonyl-", 
                     "Dichloroacetaldehyde", "Linalyl acetate", "Beta-Ocimene", 
                     "2-Thiophenecarboxylic acid", "1-Pent-3-ynylcyclopenta-1,3-diene", 
                     "1,5,6,7-Tetramethylbicyclo[3.2.0]hepta-2,6-diene", 
                     "4'-Ethylacetophenone", "Butyl citrate", 
                     "1-Methyl-4-phenyl-5-Thioxo-1,2,4-triazolidin-3-one", 
                     "Oleamide")            

all_chemicals = c("Pivalaldehyde, semicarbazone", "2-Butenal", 
                  "2,4-Dimethyl-1-heptene","4-hydroxy-4-methyl-5-phenoxy-2-Pentanone", 
                  "1-Octen-3-ol", "Butanal", 
                  "4,5-Dichloro-1,3-dioxolan-2-one", "3-Octanone", "Decane", 
                  "Mesitylene", "Benzene, 1,2,4-trimethyl-", "D-Limonene", 
                  "2-Phenylbutan-1-Ol", "1,4-diethyl-2-Methylbenzene", 
                  "1,2-diethyl-3-Methylbenzene", "1,3,8-p-Menthatriene", 
                  "4-[Dichloromethyl]-2-[[2-[1-Methyl-2-Pyrrolidinyl]Ethyl]Amino-6-Trichloromethylpyrimidine", 
                  "Benzene, (2-methyl-1-propenyl)-", "1-Phenyl-1-butene", 
                  "Linalool", "Nonanal", "2-Thiophenecarboxylic acid, 5-nonyl-", 
                  "Dichloroacetaldehyde", "Linalyl acetate", "Beta-Ocimene", 
                  "2-Thiophenecarboxylic acid", "1-Pent-3-ynylcyclopenta-1,3-diene", 
                  "1,5,6,7-Tetramethylbicyclo[3.2.0]hepta-2,6-diene", 
                  "4'-Ethylacetophenone", "Butyl citrate", 
                  "1-Methyl-4-phenyl-5-Thioxo-1,2,4-triazolidin-3-one", 
                  "Oleamide","(E)-11-hexadecenal", "(E,E)-10,12-hexadecadienal","(E)-10-hexadecenal", 
                  "(Z)-10-hexadecenal", "(E,E)-10,12-hexadecadienal", 
                  "(Z)-11-hexadecenyl acetate", "(Z)-11-hexadecenal", 
                  "(Z)-11-hexadecenol", "(4aS,7S,7aR)-nepetalactone", "(1R,4aS,7S,7aR)-nepetalactol",
                  "(E,Z)-3,13-octadecadienyl acetate", "(Z,Z)-3,13-octadecadienyl acetate", 
                  "(E)-11-tetradecenyl acetate", "(E,E)-9,11-tetradecadienyl acetate", "(E)-11-tetradecenol",
                  "(E)-11-hexadecenyl acetate", "(Z,E)-9,12-tetradecadienyl acetate", 
                  "(Z)-9-tetradecenyl acetate", "(Z)-11-hexadecenyl acetate",
                  "(Z,E)-9,12-tetradecadienol", "(Z)-9-tetradecenol", "(Z)-11-hexadecenol",
                  "(Z)-9-tetradecenyl acetate", "(Z)-9-tetradecenol", "tetradecyl acetate",
                  "(R)-delta-heptalactone", "(Z,E)-5,7-dodecadienol", 
                  "(Z,E)-5,7-dodecadienyl propionate","2,4-D","benefin","bensulide",
                  "bentazon","bromoxynil","cacodylic acid",
                  "calcium acid methanearsonate","carfentrazone","chlorsulfuron",
                  "clethodim","CHEMBL2372957","dicamba",
                  "dichlobenil","dimethenamid-P",
                  "diquat","dithiopyr","EPTC",
                  "fluazifop","fluroxypyr","foramsulfuron",
                  "glufosinate","glyphosate","halosulfuron",
                  "hexazinone","imazapyr","isoxaben",
                  "MCPA","mecoprop","metolachlor",
                  "MSMA","napropamide","oryzalin",
                  "oxadiazon","oxyfluorfen","pelargonic acid",
                  "pendimethalin","prodiamine","pronamide",
                  "quinclorac","sethoxydim","siduron",
                  "sulfosulfuron","tebuthiuron","triclopyr",
                  "trifloxysulfuron-sodium","trifluralin","abamectin", "acephate","acetamiprid",
                  "allethrin","arsenic trioxide","azadirachtin",
                  "bifenthrin","borate",
                  "carbaryl","chlorantraniliprole","clothianidin",
                  "cryolite","cyfluthrin","cypermethrin",
                  "diflubenzuron","dinotefuran","disulfoton",
                  "emamectin B1a","fipronil","fluvalinate",
                  "hydramethylnon","imidacloprid",
                  "lambda-cyhalothrin","malathion","salannin",
                  "permethrin","pyrethrin","resmethrin",
                  "rotenone","silica","spinosad",
                  "sulfluramid","thiamethoxam")

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


#########################################################################################################
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

########################################################################################################################
####                   Internal Standard - standardify(), standardifyIt()                                           ####
########################################################################################################################

standardify = function(u, v, w, x, y, z){
  standardised = (u/v*w)/x/y/z
}

################################################################################################

standardifyIt = function(data_in, int_standard, IS_ng, IS_uL, collect_time, sample_amt){
  data_in_standard = data_in[data_in$cmps_correct != int_standard, c("cmps_correct","mean_RT")]
  i = 3
  for(i in 3:ncol(data_in)){
    Input = data_in[,i][data_in$cmps_correct != int_standard]
    IS_avg = mean(data_in[data_in$cmps_correct == int_standard,,drop=T])
    IS = sum(data_in[,i][data_in$cmps_correct == int_standard])
    if(IS > 0){
      tmp = standardify(u = Input, 
                        v = IS, 
                        w = IS_ng, 
                        x = IS_uL, 
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

######################################################################################
#####################     Start of function - categorate()       #####################
######################################################################################
# i = 8
categorate = function(compounds){
  compound_cids = get_cid(compounds)
  compound_cids = compound_cids[!is.na(compound_cids$cid),]
  
  compound_SDFs = dictionate(make.unique(compounds))
  CMP_info_df = data.frame(matrix(ncol = 7, nrow = 0))
  
  Chem_data_source = c("reactives_df", "LOTUS_df",
                       "KEGG_df", "FEMA_df", 
                       "FDA_SPL_df", "Chemical")
  colnames(CMP_info_df) = Chem_data_source
  
  functional_groups = c("Isoprene", "Benzene", "Acetic acid", 
                        "Propylbenzene", "Propionic acid", 
                        "2H-furan-5-one")
  functional_identities = c("Terpenoid", "Aromatics", "FattyAcid", 
                            "Phenylpropanoid", "Polyketide", 
                            "Strigolactone", "Chemical")
  
  functional_SDFs = dictionate(functional_groups)
  functional_df = data.frame(matrix(ncol = 7, nrow = 0))
  colnames(functional_df) = functional_identities
  
  for (i in 1:length(compound_cids$cid)){
    current_cid = toString(compound_cids$cid[i])
    current_cid = gsub("[c\\\"() ]","",current_cid)
    Chemical = compound_cids$query[i]
    
    reactives_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',current_cid,'/JSON?heading=Reactive+Group')
    LOTUS_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',current_cid,'/JSON?source=LOTUS+-+the+natural+products+occurrence+database')
    KEGG_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',current_cid,'/JSON?source=KEGG')
    FEMA_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',current_cid,'/JSON?source=Flavor+and+Extract+Manufacturers+Association+(FEMA)')
    MeSH_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',current_cid,'/JSON?source=Medical+Subject+Headings+(MeSH)')
    FDA_SPL_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',current_cid,'/JSON?source=FDA/SPL+Indexing+Data')
    # CAMEO_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',current_cid,'/JSON?source=CAMEO+Chemicals')
    
    reactives_info = tryCatch(jsonlite::fromJSON(reactives_url), error = function(error) {return("None")})
    LOTUS_info = tryCatch(jsonlite::fromJSON(LOTUS_url), error = function(error) {return("None")})
    KEGG_info = tryCatch(jsonlite::fromJSON(KEGG_url), error = function(error) {return("None")})
    FEMA_info = tryCatch(jsonlite::fromJSON(FEMA_url), error = function(error) {return("None")})
    MeSH_info = tryCatch(jsonlite::fromJSON(MeSH_url), error = function(error) {return("None")})
    FDA_SPL_info = tryCatch(jsonlite::fromJSON(FDA_SPL_url), error = function(error) {return("None")})
    # CAMEO_info = tryCatch(jsonlite::fromJSON(CAMEO_url), error = function(error) {return("None")})
    
    if (reactives_info != "None"){
      reactives_value = as.matrix(unlist(reactives_info))
    } else {reactives_value = "None"}
    if (!is.na(reactives_value[2])){
      chem_name1 = paste0('"',Chemical,'"')
      chem_name1 = gsub("[c\\\"() ]", "", chem_name1)
      chem_name2 = paste0('"', toupper(Chemical), '"')
      chem_name2 = gsub("[c\\\"() ]", "", chem_name2)
      exclude_these = c("CAMEO Chemicals","Safety and Hazards",
                        "Stability and Reactivity","Reactive Group",
                        "CID", chem_name1, chem_name2)
      exclude_that = c("https:")
      
      reactives_0 = matrix(reactives_value[nchar(reactives_value[,1]) < 100])
      reactives_1 = matrix(reactives_0[!reactives_0[,1] %in% exclude_these])
      reactives_2 = matrix(reactives_1[!substr(reactives_1[,1],1,6) %in% exclude_that])
      reactives_final = c(reactives_2[!grepl("\\d",reactives_2[,1])])
    } else {reactives_final = "None"}
    reactives_df = reactives_final
    
    if (LOTUS_info != "None"){
      LOTUS_value = as.matrix(unlist(LOTUS_info))
    } else {LOTUS_value = "None"}
    if (!is.na(LOTUS_value[2])){
      LOTUS_keeps = c("Record.Reference.SourceID",
                      "Record.Reference.SourceID1",
                      "Record.Reference.SourceID2")
      LOTUS_final = LOTUS_value[row.names(LOTUS_value) %in% LOTUS_keeps]
    } else {LOTUS_final = "None"}
    LOTUS_df = LOTUS_final
    
    if (KEGG_info != "None"){
      KEGG_value = as.matrix(unlist(KEGG_info))
    } else {KEGG_value = "None"}
    if (!is.na(KEGG_value[2])){
      KEGG_final = KEGG_value[substr(KEGG_value,1,5)=="KEGG:"]
    } else {KEGG_final = "None"}
    KEGG_df = KEGG_final
    
    if (FEMA_info != "None"){
      FEMA_value = as.matrix(unlist(FEMA_info))
    } else {FEMA_value = "None"}
    if (!is.na(FEMA_value[2])){
      FEMA_final = FEMA_value[row.names(FEMA_value)=="Record.Section.Section.Information.Value.StringWithMarkup.String"]
    } else {FEMA_final = "None"}
    FEMA_df = unlist(strsplit(FEMA_final,","))
    
    if (MeSH_info != "None"){
      MeSH_value = as.matrix(unlist(MeSH_info))
    } else {MeSH_value = "None"}
    if (!is.na(MeSH_value[2])){
      MeSH_final = MeSH_value[row.names(MeSH_value)=="Record.Section.Section.Information.Value.StringWithMarkup.String"]
    } else {MeSH_final = "None"}
    MeSH_df = MeSH_final
    
    if (FDA_SPL_info != "None"){
      FDA_SPL_value = as.matrix(unlist(FDA_SPL_info))
    } else {FDA_SPL_value = "None"}
    if (!is.na(FDA_SPL_value[2])){
      FDA_SPL_final = FDA_SPL_value[row.names(FDA_SPL_value)=="Record.Reference.Name"]
    } else {FDA_SPL_final = "None"}
    FDA_SPL_df = FDA_SPL_final
    
    CMP_info_row = as.data.frame(t(rbind.data.frame(reactives_df, LOTUS_df, KEGG_df, FEMA_df, FDA_SPL_df, Chemical)))
    colnames(CMP_info_row) = Chem_data_source
    CMP_info_df = rbind(CMP_info_df, CMP_info_row)
    row.names(CMP_info_df) = NULL
  }
  
  SDF_info_df = data.frame(matrix(ncol = 9, nrow = 0))
  SDF_columns = c("Chemical", "MW", "MF", "Rings", 
                  "Groups", "GroupCounts", "Atom", 
                  "AtomCounts", "NCharges")
  
  colnames(SDF_info_df) = SDF_columns
  
  for(SDF in 1:length(compound_SDFs)){
    
    Ncharges = tryCatch(sapply(bonds(compound_SDFs[SDF], type="charge"),length), error = function(error) {return("None")})
    atom_counts = tryCatch(atomcountMA(compound_SDFs[SDF], addH = FALSE), error = function(error) {return("None")})
    group_counts = tryCatch(groups(compound_SDFs[SDF], type = "countMA"), error = function(error) {return("None")})
    ring_counts = tryCatch(rings(compound_SDFs[SDF], upper = 6, type = "count", arom = TRUE, type = "countMA"), error = function(error) {return("None")})
    Chemical = compound_SDFs[[SDF]][[4]][[9]]
    MF = MF(compound_SDFs[SDF], addH=FALSE)
    MW = MW(compound_SDFs[SDF], addH=FALSE)
    
    atom_count_atoms = colnames(atom_counts)
    group_count_groups = rownames(as.matrix(unlist(group_counts)))
    if(is.null(group_count_groups)){
      group_count_groups = "None"
    }
    SDF_info_row = as.data.frame(t(rbind.data.frame(as.vector(Chemical), 
                                                    as.vector(MW), as.vector(MF), 
                                                    as.vector(ring_counts),
                                                    group_count_groups,
                                                    as.vector(group_counts), 
                                                    atom_count_atoms,
                                                    as.vector(atom_counts), 
                                                    as.vector(Ncharges))))
    batch_test_set = fmcsBatch(compound_SDFs[SDF], functional_SDFs)
    
    if (max(as.vector(batch_test_set[,5])) > 0.95)
    {
      functional_match = names(batch_test_set[,5][batch_test_set[,5] > 0.95])
      match_SDFs = functional_SDFs[functional_SDFs@ID %in% functional_match]
      functional_match = as.integer(gsub("CMP", "", functional_match))
      functional_row = data.frame(t(rep("No", ncol(functional_df))))
      colnames(functional_row) = functional_identities
      functional_row$Chemical = Chemical
      functional_row[functional_match] = "Yes"
      functional_df = rbind(functional_df, functional_row)
      row.names(functional_df) = NULL
    }else
    {
      functional_row = data.frame(t(rep("No", ncol(functional_df))))
      colnames(functional_row) = functional_identities
      functional_row$Chemical = Chemical
      functional_df = rbind(functional_df, functional_row)
      row.names(functional_df) = NULL
    }
    
    
    colnames(SDF_info_row) = SDF_columns
    SDF_info_df = rbind(SDF_info_df, SDF_info_row)
    row.names(SDF_info_df) = NULL
  }
  data_list = list(CMP_info_df, SDF_info_df, functional_df)
  names(data_list) = c("Databases", "FMCS", "FunctionalGroups")
  return(data_list)
}


#*** PROPRIETARY - DO NOT DISTRIBUTE **