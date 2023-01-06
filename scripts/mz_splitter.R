#'[MAY BE OBSOLETE!!] - post alignment processing: mz_splitter()
#'
#'@description Reads nested list output from GCalignR row by row and 
#'splits overlapping peak identities based on a user specified number
#'of times the same m/z occurs
#'
#'@param input Nested list containing output from GCalignR
#'
#'@param samples_with User specified value for the number of times an 
#'m/z should appear on the same row to be counted as a unique compound
#'identity
#'
#'

mz_splitter = function(input, samples_with){
 floor_dec = function(x, level=1) round(x - 5*10^(-level-1), level)
 
 Orig=input$aligned$Component.Area
 # vlookup(Orig[1,3],cmp_tmp, result_column = 2, lookup_column = 1)
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
 # cmp_list2 = probs_list2 = mz_list2 = rt_list2 = list()
 rt_column = c()
 
 # i = 19
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
  # cmp_multiNames
  
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
  # probs_multiples = probs_counts$probs_current[which(probs_counts$count >= 2 & probs_counts$probs_current != 0)]
  # cmp_multiples = cmp_counts$cmp_current[which(as.numeric(as.vector(cmp_counts$count)) >= 2 & cmp_counts$cmp_current != 0)]
  # rt_multiples = rt_counts$rt_current[which(rt_counts$count >= 2 & rt_counts$rt_current != 0)]
  # area_multiples = area_counts$area_current[which(area_counts$count >= 2 & area_counts$area_current != 0)]
  
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
  # probs_solos = probs_counts$probs_current[which(probs_counts$count < 2 & probs_counts$probs_current != 0)]
  # cmp_solos = cmp_counts$cmp_current[which(as.numeric(as.vector(cmp_counts$count)) < 2 & cmp_counts$cmp_current != 0)]
  # rt_solos = rt_counts$rt_current[which(rt_counts$count < 2 & as.numeric(as.vector(rt_counts$rt_current)) != 0)]
  # area_solos = area_counts$area_current[which(area_counts$count < 2 & as.numeric(as.vector(area_counts$area_current)) != 0)]
  
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
                                   
                                   max_location = solos_matrix$mz %in% mz_solos
                                   max_prob = max(as.numeric(as.vector(unlist(solos_matrix$Prob[max_location]))))
                                   probs_row[probs_no] <- 0
                                   probs_row[probs_yes] <- max_prob  
                                   
                                   best_cmp = solos_matrix$Compound[solos_matrix$Prob == max_prob]
                                   cmp_row[cmp_no] <- "NA"
                                   cmp_row[cmp_yes] <- as.character(unlist(best_cmp))
                                   
                                   best_mz = solos_matrix$mz[solos_matrix$Prob == max_prob]
                                   mz_row[mz_no] <- 0
                                   mz_row[mz_yes] <- as.numeric(unlist(best_mz))
                                   
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
                                   # mz_row = mz_matrix[i,]
                                   # probs_row = probs_matrix[i,]
                                   # cmp_row = cmp_matrix[i,]
                                   # rt_row = rt_matrix[i,]
                                   # area_row = area_matrix[i,]
                                   # 
                                   # cmp_matrix_none = rbind(cmp_matrix_none, 
                                   #                         cmp_row)
                                   # probs_matrix_none = rbind(probs_matrix_none, 
                                   #                           probs_row)
                                   # mz_matrix_none = rbind(mz_matrix_none, 
                                   #                        mz_row)
                                   # rt_matrix_none = rbind(rt_matrix_none, 
                                   #                        rt_row)
                                   # area_matrix_none = rbind(area_matrix_none, 
                                   #                          area_row)
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