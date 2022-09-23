#'Intelligent row merging function that checks within a range of likely 
#'RT overlaps. 
#'
#'@description Likely RT overlaps are calculated from the input data. 
#'Instances where multiple compounds are tentatively identified 
#'within that range are compared in terms of: a) number of instances a 
#'compound was identified across all samples; b) maximum match probability 
#'for each uniquely identified chemical; and, c) molecular mass. The final 
#'weighting metric involves a combination of a, b, and c. The best compound 
#'identification is determined by this score, component areas across the RT 
#'range are summed and the lowest RT is selected. A few more cleaning steps 
#'are required to remove duplicates and missing values then the merged data 
#'set is returned.
#'
#'@param input_list A list containing our current "aligned_peaks" object 
#'equivalent.
#'

theMerger = function(input_list){
 probs_matrix = input_list$MatchFactor
 cmp_matrix = input_list$Compounds
 RT_matrix = input_list$RT
 mz_matrix = input_list$MZ
 
 best_cmps = vector()
 probs_cmps = vector()
 mean_RT = vector()
 filler_df = data.frame(matrix(ncol = 3, 
                               nrow = 0))
 colnames(filler_df) = c("Probability",
                         "Compound",
                         "RT")
 
 for (i in 1:nrow(probs_matrix)){
  probs_row = probs_matrix[i,][!is.na(probs_matrix[i,])]
  probs_row = probs_row[probs_row != "0"]
  row.names(probs_row) = NULL
  colnames(probs_row) = NULL
  if (length(probs_row) > 1){
   probs_row = max(probs_row)
  } else {}
  
  cmp_row = cmp_matrix[i,][!is.na(cmp_matrix[i,])]
  cmp_row = cmp_row[cmp_row != "0"]
  row.names(cmp_row) = NULL
  colnames(cmp_row) = NULL
  if (length(cmp_row) > 1){
   prob_again = probs_matrix[i,][!is.na(probs_matrix[i,])]
   prob_again = prob_again[prob_again != "0"]
   cmp_picker = matrix(ncol = 2, 
                       nrow = length(cmp_row))
   cmp_picker[,1] = prob_again
   cmp_picker[,2] = cmp_row
   colnames(cmp_picker) = c("prob", "cmp")
   cmp_row = cmp_picker[,2][cmp_picker[,1] == max(cmp_picker[,1])]
  }
  
  RT_row = RT_matrix[i,][!is.na(RT_matrix[i,])]
  RT_row = RT_row[RT_row != "0"]
  row.names(RT_row) = NULL
  colnames(RT_row) = NULL
  if (length(RT_row) > 1){
   RT_row = min(RT_row)
  } else {}
  filler_df[i,] = c(as.numeric(max(probs_row)), 
                    as.character(cmp_row), 
                    as.numeric(min(RT_row)))
 }
 aligned_with_CMPS = as.data.frame(cbind(filler_df, 
                                         input_list$Area))
 unique_CMPs = unique(aligned_with_CMPS$Compound)
 num_unique_CMPs = length(unique(aligned_with_CMPS$Compound))
 CMP_counts = c()
 CMP_mass = c()
 CMP_name_SDF = c()
 all_dat = c()
 
 print("Removing any compounds that only occur in 1 sample.")
 for (single in 1:num_unique_CMPs){
   current_CMP = unique_CMPs[single]
   CMP_count = nrow(aligned_with_CMPS[aligned_with_CMPS$Compound == current_CMP,])
   CMP_counts = c(CMP_counts, 
                  CMP_count)
 }
 CMP_count_df = data.frame(cbind(unique_CMPs, CMP_counts))
 single_matches = paste0(CMP_count_df$unique_CMPs[CMP_count_df$CMP_counts == 1])
 aligned_no_solos = aligned_with_CMPS[!aligned_with_CMPS$Compound %in% single_matches,]
 
 unique_CMPs = unique(aligned_no_solos$Compound)
 num_unique_CMPs = length(unique(aligned_no_solos$Compound))
 CMP_counts = c()
 CMP_mass = c()
 CMP_name_SDF = c()
 all_dat = c()
 
 # cmp = 41
 print("Acquiring exact mass data for chemicals published on PubChem.")
 for (cmp in 1:num_unique_CMPs){
  current_CMP = unique_CMPs[cmp]
  CMP_count = nrow(aligned_no_solos[aligned_no_solos$Compound == current_CMP,])
  CMP_counts = c(CMP_counts, 
                 CMP_count)
  current_cid = get_cid(current_CMP)
  current_cid = toString(current_cid[[1,2]])
  current_cid = gsub("NA", 
                     "180", 
                     current_cid)
  current_cid = gsub("[c\\\"() ]",
                     "",
                     current_cid)
  current_cid = gsub("\n", 
                     "", 
                     current_cid)
  mass_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
                    current_cid,
                    '/SDF')
  mass_info = tryCatch(read.SDFset(mass_url), 
                       error = function(error) {return("None")})
  mass_tmp = mass_info@SDF[[1]]@datablock[names(mass_info@SDF[[1]]@datablock) == "PUBCHEM_EXACT_MASS"]
  CMP_name_tmp = mass_info@SDF[[1]]@datablock[names(mass_info@SDF[[1]]@datablock) == "PUBCHEM_IUPAC_CAS_NAME"]
  if (length(CMP_name_tmp) <1){
   CMP_name_tmp = current_CMP
  } else {}
  
  names(mass_tmp) = NULL
  names(CMP_name_tmp) = NULL
  all_tmp = paste0(CMP_count, 
                   ';', 
                   mass_tmp, 
                   ';', 
                   CMP_name_tmp, 
                   ';', 
                   current_CMP)
  print(paste0(cmp, ':', all_tmp))
  all_dat = c(all_dat, 
              all_tmp)
 }
 print("Performing the first/easiest merge (identical chemical names).")
 dat_list = strsplit(all_dat, 
                     ";")
 dat_df = do.call(rbind.data.frame, 
                  dat_list)
 colnames(dat_df) = c("Count", 
                      "Mass", 
                      "Chemical1", 
                      "Chemical2")
 area_aggregate = aligned_no_solos[,-c(1,3)]
 area_aggregate[is.na(area_aggregate)] = 0
 cols_num = c(2:ncol(area_aggregate))
 area_aggregate[cols_num] = sapply(area_aggregate[cols_num], as.numeric)
 area_dat = aggregate(. ~ Compound, 
                      data = area_aggregate, 
                      FUN = sum)
 RT_dat = aggregate(RT ~ Compound, 
                    data = aligned_no_solos[,c(2,3)], 
                    FUN = min)
 prob_dat = aggregate(Probability ~ Compound, 
                      data = aligned_no_solos[,c(2,1)], 
                      FUN = max)
 RT = as.numeric(RT_dat[,-1])
 prob = as.numeric(prob_dat[,-1])
 
 
 mass_column = rep("NA", 
                   nrow(area_dat))
 count_column = rep("NA", 
                    nrow(area_dat))
 
 # chems = 1867
 for (chems in 1:nrow(area_dat)){
  match_chem1 = as.character(dat_df$Chemical1[chems])
  match_chem2 = as.character(dat_df$Chemical2[chems])
  match_chems = c(match_chem1, 
                  match_chem2)
  logical_chem1 = match_chem1 %in% area_dat$Compound
  logical_chem2 = match_chem2 %in% area_dat$Compound
  
  if (logical_chem1 | logical_chem2){
   mass_matched = as.numeric(paste0(dat_df$Mass[chems]))
   count_matched = as.numeric(paste0(dat_df$Count[chems]))
   mass_column[area_dat$Compound %in% match_chems] = mass_matched
   count_column[area_dat$Compound %in% match_chems] = count_matched
  } else {}
 }
 
 first_merger = data.frame(cbind(RT, 
                                 prob, 
                                 count_column, 
                                 mass_column, 
                                 area_dat))
 first_merger = first_merger[order(first_merger$RT),]
 
 RT_start = min(as.numeric(aligned_no_solos$RT))
 RT_end = max(as.numeric(aligned_no_solos$RT))
 
 RT_range = (RT_end-RT_start)/num_unique_CMPs
 
 all_df = data.frame(matrix(ncol = ncol(area_dat[,-1])+3, 
                            nrow = 0), 
                     stringsAsFactors = FALSE)
 # all_df[is.na(all_df)] = 0
 colnames(all_df) = c("RT", 
                      "Mass", 
                      "Chemical", 
                      colnames(area_dat[,-1]))
 
 merger_rows = first_merger
 print("Revving up the heavier merge that runs on science.")
 for (RT_step1 in 1:nrow(first_merger)){
  RT_1 = merger_rows$RT[RT_step1]
  RT_2 = RT_1 + RT_range
  RT_rows = merger_rows$RT >= RT_1 & merger_rows$RT <= RT_2
  
  df_filler = data.frame(merger_rows[RT_rows,], stringsAsFactors = FALSE)
  
  if (nrow(df_filler) > 1){
   score = c()
   for (i in 1:nrow(df_filler)){
    prob_tmp = (-1+2.7183^(df_filler$prob[i]/100))
    count_tmp = 2-1/as.numeric(paste0(df_filler$count_column[i]))
    score_tmp = (prob_tmp/as.numeric(paste0(df_filler$count_column[i])))*count_tmp
    score = c(score, 
              score_tmp)
   }
   df_RT_range = data.frame(cbind(score, 
                                  df_filler), 
                            stringsAsFactors = FALSE)
   mass_matrix = data.frame(matrix(nrow = nrow(df_filler), 
                                   ncol = nrow(df_filler)), 
                            stringsAsFactors = FALSE)
   colnames(mass_matrix) = paste0("CMP", 
                                  rep.int(1:nrow(df_filler),
                                          1))
   rownames(mass_matrix) = paste0("CMP", 
                                  rep.int(1:nrow(df_filler),
                                          1))
   for (m in 1:nrow(df_filler)){
    for (n in 1:nrow(df_filler)){
     diffs = all.equal(as.numeric(paste0(df_RT_range$mass_column[m])), 
                       as.numeric(paste0(df_RT_range$mass_column[n])))
     diffs = tryCatch(strsplit(diffs, ": "), 
                      error = function(error) {return(NA)})
     diffs = suppressWarnings(as.numeric(paste0(diffs[[1]][2])))
     mass_matrix[m,n] = diffs
    }
   }
   mass_diff_weights = as.numeric(paste0(colSums(mass_matrix, 
                                                 na.rm = T)))
   df_RT_range = data.frame(cbind(mass_diff_weights, 
                                  df_RT_range), 
                            stringsAsFactors = FALSE)
   combo_score = df_RT_range$score/df_RT_range$mass_diff_weights
   df_RT_range = data.frame(cbind(combo_score, 
                                  df_RT_range), 
                            stringsAsFactors = FALSE)
   area_summed = as.numeric(paste0(colSums(df_RT_range[,-c(1:8)])))
   # if(IS %in% df_RT_range$Compound){
   #   best_CMP = as.character(paste0(IS))
   #   # if(paste0(all_df$Chemical[nrow(all_df)], "test") == paste0(best_CMP, "test")) next
   # } else {
     best_CMP = as.character(paste0(df_RT_range$Compound[df_RT_range$combo_score == max(df_RT_range$combo_score)]))
   # }
   if(paste0(all_df$Chemical[nrow(all_df)], "test") == paste0(best_CMP, "test")) next
   best_RT = as.numeric(paste0(min(df_RT_range$RT)))
   best_mass = as.numeric(paste0(df_RT_range$mass_column[df_RT_range$combo_score == max(df_RT_range$combo_score)]))
   
   all_tmp = cbind(best_RT, 
                   best_mass, 
                   best_CMP, 
                   rbind(area_summed))
   row.names(all_tmp) = NULL
   colnames(all_tmp) = c("RT", 
                         "Mass", 
                         "Chemical", 
                         colnames(df_RT_range[,-c(1:8)]))
  } else {
   row_tmp = data.frame(first_merger[RT_rows,], 
                        stringsAsFactors = FALSE)
   row.names(row_tmp) = NULL
   all_tmp = row_tmp[,-c(2,3)]
   colnames(all_tmp) = c("RT", 
                         "Mass", 
                         "Chemical", 
                         colnames(df_RT_range[,-c(1:8)]))
  }
  if(length(all_df$Chemical[nrow(all_df)]) < 1 & RT_step1 != 1) next
  all_df = rbind(all_df, 
                 all_tmp)
 }
 all_df_final = all_df[!is.na(all_df$RT),]
 unduped_dat = all_df_final[,-c(1:3)]
 print("Final tidying of merged data. Almost ready.")
 for (col_final in 1:ncol(unduped_dat)){
  unduped_dat[,col_final][duplicated(unduped_dat[,col_final]) & unduped_dat[,col_final] != 0] = 0
 }
 all_df_last_merge = data.frame(cbind(all_df_final[,c(1:3)],unduped_dat))
 
 if(any(duplicated(all_df_last_merge$Mass))){
  all_df_last_merge$Mass[duplicated(all_df_last_merge$Mass)] = NA
 }
 all_df_last_merge = all_df_last_merge[!is.na(all_df_last_merge$Mass),]
 
 empty_rows = c()
 for(row in 1:nrow(all_df_last_merge)){
  test_row = as.numeric(paste0(unlist(all_df_last_merge[row,-c(1:3)])))
  test_sum = sum(test_row)
  empty_rows = c(empty_rows, test_sum)
 }
 all_df_last_merge = data.frame(cbind(empty_rows, 
                                      all_df_last_merge))
 all_df_true_final = all_df_last_merge[all_df_last_merge$empty_rows != 0,-1]
 return(all_df_true_final)
}
