#'Post Alignment Toolbox
#'
#'@description Function that allows the user of the package to take
#'control of their data by subsetting according to user defined 
#'parameters. 
#'
#'Decontaminate, a function to eliminate compounds that are not of
#'interest from your dataset based on a user-defined character string option 
#'for exact or non-exact matches, default is exact.

data_in = unknowns_spread
# chemicals = new_chems
RT_cutoff = 17
RT_search_range = 2
# IS = "Tetradecane"
chemicals = c("Acetic ester", "Ethyl hexanoate", "Octanal", "Undecane", "Methyl salicylate")

mzExacto <- function(data_in, chemicals, RT_cutoff = 0, RT_search_range = 2){
  area_matrix = data_in$Area
  probs_matrix = data_in$MatchFactor
  cmp_matrix = data_in$Compounds
  RT_matrix = data_in$RT
  mz_matrix = data_in$MZ
  mass_matrix = data_in$Mass
  rtBYmass_matrix = data_in$rtBYmass
  
  input_CMPs_long = c()
  
  # this needs to be for the list of chemicals
  # w = NULL
  
  for (w in 1:nrow(cmp_matrix)){
    input_CMPs_tmp = cmp_matrix[w,][!is.na(cmp_matrix[w,])]
    input_CMPs_long = c(input_CMPs_long, input_CMPs_tmp)
  }
  mz_inputs_tmp = as.data.frame(mz_matrix)
  # row.names(mz_inputs_tmp) = input_CMPs_long
  mz_inputs = mz_inputs_tmp[input_CMPs_long %in% chemicals,]
  input_MZs_long = c()
  
  for (x in 1:nrow(mz_inputs)){
    input_MZs_tmp = mz_inputs[x,][!is.na(mz_inputs[x,])]
    input_MZs_long = c(input_MZs_long, input_MZs_tmp)
  }
  
  # RTs_tmp = RT_matrix[z,][!is.na(RT_matrix[z,])]
  # RTs_long = c(RTs_long, RTs_tmp)
  # 
  # masses_tmp = mass_matrix[z,][!is.na(mass_matrix[z,])]
  # masses_long = c(masses_long, masses_tmp)
  # 
  # MZs_tmp = mz_matrix[z,][!is.na(mz_matrix[z,])]
  # MZs_long = c(MZs_long, MZs_tmp)
  # 
  # rtBYmass_tmp = rtBYmass_matrix[z,][!is.na(rtBYmass_matrix[z,])]
  # rtBYmass_long = c(rtBYmass_long, rtBYmass_tmp)
  # ###############################################################
  
  
  df_columns = ncol(area_matrix)+3
  df_rows = length(chemicals)
  all_df_last_merge = data.frame(matrix(nrow = df_rows, ncol = df_columns))
  num_unique_CMPs = length(chemicals)
  
  final_list = list()
  
  colnames(all_df_last_merge) = c("RT", "Mass", "Chemical", colnames(area_matrix))
  
  chem_mass_set = c()
  mz_mass_list = list()
  
  cat("Acquiring relevant info from Pubchem [pubchem.ncbi.nlm.nih.gov]. Please be patient!!")
  for (chem in 1:length(chemicals)){
    
    alt_trigger = F
    current_CMP = chemicals[chem]
    chem_cid = get_cid(chemicals[chem])
    chem_cid = paste0(chem_cid[[1,2]])
    
    chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid,'/JSON?heading=GC-MS')
    
    chem_info = tryCatch(jsonlite::fromJSON(chem_url), 
                         error = function(error) {return("None")})
    if (chem_info == "None"){
      chem_cid = 180
      
      chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid,'/JSON?heading=GC-MS')
      
      chem_info = tryCatch(jsonlite::fromJSON(chem_url), 
                           error = function(error) {return("None")})
      # mz_primary_alt = "1a"
      alt_trigger = T
      } else{}
    
    mass_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
                      chem_cid,
                      '/SDF')
    mass_info = tryCatch(read.SDFset(mass_url), 
                         error = function(error) {return("None")})
    mass_tmp = mass_info@SDF[[1]]@datablock[names(mass_info@SDF[[1]]@datablock) == "PUBCHEM_EXACT_MASS"]
    CMP_name_tmp = mass_info@SDF[[1]]@datablock[names(mass_info@SDF[[1]]@datablock) == "PUBCHEM_IUPAC_CAS_NAME"]
    if (length(CMP_name_tmp) <1){
      CMP_name_tmp = current_CMP
    } else {}
    
    chem_json_dat = data.frame(unlist(chem_info$Record$Section))
    colnames(chem_json_dat) = "info"
    mz_rows_1 = grep("[[:digit:]]+\\.[[:digit:]]\\ [[:digit:]]+\\.[[:digit:]]", chem_json_dat$info)
    mz_rows_2 = grep("[[:digit:]]{2,3}\\ [[:digit:]]+\\.[[:digit:]]", chem_json_dat$info)
    mz_rows_3 = grep("[[:digit:]]{2,3}\\ [[:digit:]]{2,3}", chem_json_dat$info)
    mz_rows_4 = grep("^[[:digit:]]{2,3}$", chem_json_dat$info)
    
    if (length(mz_rows_1) > 0) {
      mz_rows = mz_rows_1
    } else {
      if (length(mz_rows_2) > 0) {
        mz_rows = mz_rows_2
      } else {
        if(length(mz_rows_3) > 0) {
          mz_rows = mz_rows_3
        } else {
          if(length(mz_rows_4) > 0) {
          mz_rows = mz_rows_4
        } else {
          mz_rows = 0
          alt_trigger = T
        }
        }
      }
    }
    
    
    chem_mz_matches = unique(paste0(unlist(strsplit(paste0(chem_json_dat[mz_rows,])," "))))
    chem_mz_matches2 = gsub("\\.0\\>","",chem_mz_matches)
    chem_mz_matches2 = gsub("\\.[[:digit:]]0\\>","",chem_mz_matches2)
    
    mz_count = 1:length(chem_mz_matches2)
    odds = mz_count[mz_count%%2 != 0]
    evens = mz_count[mz_count%%2 == 0]
    
    mz_primary = chem_mz_matches2[odds]
    mz_secondary = chem_mz_matches2[evens]
    chem_mz_matches2 = mz_primary[order(mz_secondary, decreasing = T)]
    
    if(alt_trigger == T){
      mz_primary = input_MZs_long[input_CMPs_long == input_CMPs_tmp]
      mass_tmp = "NA"}
    
    mz_print = paste(chem_mz_matches2, collapse = " | ")
    total_chems = length(chemicals)
    
    step_printer = c(1:num_unique_CMPs)
    step_cmper = chem/25
    
    if(step_cmper %in% step_printer | chem == 1){
      cat(paste0('\n', '[Current/Total]', ' |--Exact Mass--|--CMP Name 1--|--CMP Name 2--|--Top MZ Peaks--|', '\n'))
    } else {}
    
    names(mass_tmp) = NULL
    names(CMP_name_tmp) = NULL
    all_tmp = paste0(mass_tmp, 
                     ';', 
                     CMP_name_tmp, 
                     ';', 
                     current_CMP)
    all_print = paste0(mass_tmp, 
                       ' | ', 
                       CMP_name_tmp, 
                       ' | ',
                       current_CMP,
                       ' | ',
                       mz_print)
    cat(paste0('[', chem, ' / ', num_unique_CMPs, ']', ': ', all_print, '\n'))
    
    mz_mass_list[[chem]] = list(mass_tmp, mz_primary)
  }
  if(length(mz_mass_list) < length(chemicals)){mz_mass_list[[length(chemicals)]] = "NA"}
  names(mz_mass_list) = chemicals
  
  mass_each_chem = sapply(mz_mass_list, '[[', 1)
  mz_each_chem = sapply(mz_mass_list, '[[', -1)
  mass_each_chem[mass_each_chem == "NA"] = 999.999
  
  mass_orderer = as.numeric(paste0(mass_each_chem))

  masses_ordered = mass_each_chem[order(mass_orderer)]
  mz_mass_ordered = mz_each_chem[order(mass_orderer)]
  chems_ordered = names(masses_ordered)
  
  CMPs_long = c()
  probs_long = c()
  RTs_long = c()
  masses_long = c()
  MZs_long = c()
  rtBYmass_long = c()
  MZs_published = list()
  
  # z = NULL
  for (z in 1:nrow(cmp_matrix)){
    CMPs_tmp = cmp_matrix[z,][!is.na(cmp_matrix[z,])]
    CMPs_long = c(CMPs_long, CMPs_tmp)
    
    probs_tmp = probs_matrix[z,][!is.na(probs_matrix[z,])]
    probs_long = c(probs_long, probs_tmp)
    
    RTs_tmp = RT_matrix[z,][!is.na(RT_matrix[z,])]
    RTs_long = c(RTs_long, RTs_tmp)
    
    masses_tmp = mass_matrix[z,][!is.na(mass_matrix[z,])]
    masses_long = c(masses_long, masses_tmp)
    
    MZs_tmp = mz_matrix[z,][!is.na(mz_matrix[z,])]
    MZs_long = c(MZs_long, MZs_tmp)
    
    rtBYmass_tmp = rtBYmass_matrix[z,][!is.na(rtBYmass_matrix[z,])]
    rtBYmass_long = c(rtBYmass_long, rtBYmass_tmp)
    
    chem_cid = get_cid(CMPs_tmp)
    chem_cid = paste0(chem_cid[[1,2]])
    chem_cid = gsub("NA",
                    "180",
                    chem_cid)
    
    chem_SDF_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
                    chem_cid,
                      '/SDF')
    
    chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid,'/JSON?heading=GC-MS')
    
    chem_info = tryCatch(jsonlite::fromJSON(chem_url),
                         error = function(error) {return("None")})
    
    if(chem_info == "None"){next}
    chem_json_dat = data.frame(unlist(chem_info$Record$Section))
    colnames(chem_json_dat) = "info"
    mz_rows_1 = grep("[[:digit:]]+\\.[[:digit:]]\\ [[:digit:]]+\\.[[:digit:]]", chem_json_dat$info)
    mz_rows_2 = grep("[[:digit:]]{2,3}\\ [[:digit:]]+\\.[[:digit:]]", chem_json_dat$info)
    mz_rows_3 = grep("[[:digit:]]{2,3}\\ [[:digit:]]{2,3}", chem_json_dat$info)
    mz_rows_4 = grep("^[[:digit:]]{2,3}$", chem_json_dat$info)
    # chem_json_dat$info[startsWith(chem_json_dat$info, "[0-9]{2,3}")]
    
    if (length(mz_rows_1) > 0) {
      mz_rows = mz_rows_1
    } else {
      if (length(mz_rows_2) > 0) {
        mz_rows = mz_rows_2
      } else {
        if(length(mz_rows_3) > 0) {
          mz_rows = mz_rows_3
        } else {
          if(length(mz_rows_4) >0) {
            # mz_rows_4 = mz_rows_4[-1]
            mz_rows = mz_rows_4
          }else{next}
        }
      }
    }
    #######################################################################################
    chem_mz_matches = unique(paste0(unlist(strsplit(paste0(chem_json_dat[mz_rows,])," "))))
    chem_mz_matches2 = gsub("\\.0\\>","",chem_mz_matches)
    chem_mz_matches2 = gsub("\\.[[:digit:]]0\\>","",chem_mz_matches2)
    
    mz_count = 1:length(chem_mz_matches2)
    odds = mz_count[mz_count%%2 != 0]
    evens = mz_count[mz_count%%2 == 0]
    
    mz_primary = chem_mz_matches2[odds]
    mz_secondary = chem_mz_matches2[evens]
    chem_mz_matches2 = mz_primary[order(mz_secondary, decreasing = T)]
    
    if(alt_trigger == T){
      # MZs_long[CMPs_long == CMPs_tmp]
      mz_primary = MZs_long[CMPs_long == CMPs_tmp]
      chem_mz_matches2 = mz_primary
    }
    if(chem_cid == 180){mz_primary = MZs_tmp}
    #######################################################################################
    mz_print = paste(chem_mz_matches2, collapse = " | ")
    
    step_printer = c(1:num_unique_CMPs)
    step_cmper = z/10
    
    if(step_cmper %in% step_printer | z == 1){
      cat(paste0('\n', '[Current/Total]', ' |--Chemical Name--|--Top MZ Peaks--|', '\n'))
    } else {}
    
    names(mass_tmp) = NULL
    
    all_print = paste0(CMPs_tmp,
                       ': ',
                       mz_print)
    cat(paste0('[', z, '/', nrow(cmp_matrix), ']', '-', all_print, '\n'))
    
    MZs_published[[z]] = unique(mz_primary)
  }
  
  
  if(length(MZs_published) < length(masses_long)){MZs_published[[length(masses_long)]] = NA}
  names(MZs_published) = paste0(RTs_long, ' | ', masses_long)
  MZs_probs = MZs_published
  names(MZs_probs) = paste0(RTs_long, ' | ', probs_long)
  
  all_best_RTs = c()
  previous_areas = c()
  
  row_find_df = data.frame(cbind(CMPs_long, as.numeric(paste0(RTs_long)), as.numeric(paste0(probs_long)), as.numeric(paste0(masses_long)), MZs_long))
  rownames(row_find_df) = NULL
  colnames(row_find_df) = c("Chemical", "RT", "Probs", "Mass", "MZ")
  
  probs_found = c()
  RT_found = c()
  area_found = c()
  previous_RTs = c()
  
  masses_ordered = mass_each_chem[order(mass_orderer)]
  mz_mass_ordered = mz_each_chem[order(mass_orderer)]
  chems_ordered = names(masses_ordered)
  # k = NULL
  # k = 12
  # RT_cutoff = 17
  # RT_search_range = 1
  
  # for(j in 1:length(chems_ordered)){
  #   
  # }
  RT_search_range = 0
  k = 2
  for(k in 1:length(chems_ordered)){
    row_codes = c()
    if(k > 1){previous_row = final_list[[k-1]]}else{previous_row = NA}
    # MZs_published[names(MZs_published) %in% names(MZ_possibles)]
    
    # best_RTs = lapply(seq_along(RT_df_tmp_apply),
    #                   function(i) {min(RT_df_tmp_apply[[i]][probs_df_tmp_apply[[i]] >= max(probs_df_tmp_apply[[i]], na.rm = T) - 10])})
    
    all_masses = as.numeric(paste0(masses_long))
    current_mass_tmp = as.numeric(paste0(masses_ordered[k]))
    current_mass = current_mass_tmp + (current_mass_tmp/2)
    MZ_possibles_1 = MZs_published[all_masses <= current_mass]
    MZ_probs_1 = MZs_probs[all_masses <= current_mass]
    
    MZ_possibles_split = strsplit(names(MZ_possibles_1), ' | ')
    MZ_possibles_RTs = sapply(MZ_possibles_split, "[[", 1)
    MZ_possibles_masses = sapply(MZ_possibles_split, "[[", 3)
    
    MZ_probs_split = strsplit(names(MZ_probs_1), ' | ')
    MZ_possibles_probs = sapply(MZ_probs_split, "[[", 3)
    # possibles_RTs = RTs_long[all_masses <= current_mass]
    
    if(RT_cutoff > 0){MZ_possibles = MZ_possibles_1[as.numeric(paste0(MZ_possibles_RTs )) <= RT_cutoff]}else{MZ_possibles = MZ_possibles_1}
    
    # if(k == 1){
    #   all_masses = as.numeric(paste0(masses_long))
    #   current_mass_tmp = as.numeric(paste0(masses_ordered[k]))
    #   current_mass = current_mass_tmp + (current_mass_tmp/2)
    #   MZ_possibles_1 = MZs_published[all_masses <= current_mass]
    #   
    #   MZ_possibles_bad_1 = strsplit(names(MZ_possibles_1), ' | ')
    #   MZ_possibles_bad_RT = sapply(MZ_possibles_bad_1, "[[", 1)
    #   MZ_possibles_bad_mass = sapply(MZ_possibles_bad_1, "[[", 3)
    #   MZ_possibles_min_RT_1 = sapply(MZ_possibles_bad_1, "[[", 1)
    #   
    #   # all_mass_current[p] == 58.041864811 & !between(as.numeric(paste0(all_RT_current[p])), as.numeric(paste0(all_RT_current[p])) - RT_ranger, as.numeric(paste0(all_RT_current[p])) + RT_ranger)
    #   
    #   if(RT_cutoff > 0){MZ_possibles = MZ_possibles_1[as.numeric(paste0(MZ_possibles_bad_RT)) <= RT_cutoff]}
    #   MZ_possibles_min_RT_1 = min(as.numeric(paste0(MZ_possibles_min_RT_1)))
    #   MZ_possibles = MZ_possibles[as.numeric(paste0(MZ_possibles_bad_RT)) <= as.numeric(paste0(MZ_possibles_min_RT_1)) + RT_search_range]
    #   
    # } else {
    #   all_masses = as.numeric(paste0(masses_long))
    #   previous_mass = as.numeric(paste0(masses_ordered[k-1]))
    #   if(previous_mass == 999.999){previous_mass = max(as.numeric(paste0(masses_ordered[masses_ordered < 999.99])))}
    #   # if(previous_mass == 999.999){previous_mass = as.numeric(paste0(masses_ordered[k-3]))}
    #   # if(previous_mass == 999.999){previous_mass = as.numeric(paste0(masses_ordered[k-4]))}
    #   # if(previous_mass == 999.999){previous_mass = as.numeric(paste0(masses_ordered[k-5]))}
    #   # if(previous_mass == 999.999){previous_mass = as.numeric(paste0(masses_ordered[k-6]))}
    #   # if(previous_mass == 999.999){previous_mass = as.numeric(paste0(masses_ordered[k-7]))}
    #   # if(previous_mass == 999.999){previous_mass = as.numeric(paste0(masses_ordered[k-8]))}
    #   # if(previous_mass == 999.999){previous_mass = as.numeric(paste0(masses_ordered[k-9]))}
    #   current_mass_tmp = as.numeric(paste0(masses_ordered[k]))
    #   current_mass = current_mass_tmp + (current_mass_tmp/2)
    #   
    #   
    #   MZ_possibles_1 = MZs_published[all_masses <= current_mass]
    #   
    #   MZ_possibles_bad_1 = strsplit(names(MZ_possibles_1), ' | ')
    #   MZ_possibles_bad_RT = sapply(MZ_possibles_bad_1, "[[", 1)
    #   MZ_possibles_bad_mass = sapply(MZ_possibles_bad_1, "[[", 3)
    #   MZ_possibles_min_RT_1 = sapply(MZ_possibles_bad_1, "[[", 1)
    #   RT_nums = as.numeric(paste0(MZ_possibles_bad_RT))
    #   RT_all_tmp = RTs_long[as.numeric(paste0(RTs_long)) %in% RT_nums]
    #   probs_all_tmp = probs_long[as.numeric(paste0(RTs_long)) %in% RT_nums]
    #   best_RTs_1 = RT_all_tmp[as.numeric(paste0(probs_all_tmp)) >= max(as.numeric(paste0(probs_all_tmp))) - 10]
    #   best_RT = min(as.numeric(paste0(best_RTs_1)))
    #   # best_RTs = lapply(seq_along(RT_df_tmp_apply),
    #   #                   function(i) {RT_df_tmp_apply[[i]][probs_df_tmp_apply[[i]] == max(probs_df_tmp_apply[[i]], na.rm = T)]})
    #   
    #   if(RT_cutoff > 0){MZ_possibles = MZ_possibles_1[as.numeric(paste0(MZ_possibles_bad_RT)) <= RT_cutoff]}
    #   MZ_possibles_min_RT_1 = as.numeric(paste0(MZ_possibles_min_RT_1[as.numeric(paste0(MZ_possibles_min_RT_1)) == best_RT]))
    #   
    #   if(k == 1){
    #     MZ_possibles = MZ_possibles[as.numeric(paste0(MZ_possibles_bad_RT)) <= as.numeric(paste0(best_RT)) + RT_search_range]
    #   }else{
    #     RT_previous = all_best_RTs[k-1]
    #     MZ_possibles = MZ_possibles[which(as.numeric(paste0(MZ_possibles_bad_RT)) > RT_previous & as.numeric(paste0(MZ_possibles_bad_RT)) < as.numeric(paste0(best_RT)) + RT_search_range)]
    #   }
    # }
    
    # RT_mass_list = strsplit(names(MZ_possibles),' | ')
    # 
    # 
    # all_RT_current = sapply(RT_mass_list, "[[", 1)
    # all_mass_current = sapply(RT_mass_list, "[[", 3)
    
    # RT_ranger = sd(as.numeric(paste0(all_RT_current)))/2
    # bad_mass = as.numeric(paste0(bad_mass_list[[1]][3]))
    # bad_RT = as.numeric(paste0(bad_mass_list[[1]][1]))
    
    
    # bad_mass_list = strsplit(names(MZ_possibles[p]),' | ')
    # bad_mass = as.numeric(paste0(bad_mass_list[[1]][3]))
    # bad_RT = as.numeric(paste0(bad_mass_list[[1]][1]))
    # RT_cutoff = 15
    # p = 3
    
    for(p in 1:length(MZ_possibles)){
      tryCatch(length(MZ_possibles[[p]]), error = function(error) {next})
      if(length(MZ_possibles[[p]]) == 0){next}
      
      # bad_RT = as.numeric(paste0(bad_mass_list[[1]][1]))
      # if(all_mass_current[p] == 58.041864811 & !between(as.numeric(paste0(all_RT_current[p])), as.numeric(paste0(all_RT_current[p])) - RT_ranger, as.numeric(paste0(all_RT_current[p])) + RT_ranger)){next}
      # if(as.numeric(paste0(all_RT_current[p])) >= as.numeric(paste0(RT_cutoff))){next}
      # if(bad_RT < current_RT - 3 | bad_RT > current_RT + 3){next}
      if(any(mz_mass_ordered[[k]] %in% MZ_possibles[[p]])){
        row_codes_tmp = names(MZ_possibles[MZ_possibles[[p]] %in% mz_mass_ordered[[k]]])
        row_codes = c(row_codes, row_codes_tmp)
      } else {next}
    }
    
    row_codes = row_codes[!is.na(row_codes)]
    if (is.null(row_codes)){
      row_codes_tmp = MZs_long[CMPs_long == chems_ordered[k]]
      row_codes = MZs_long == row_codes_tmp
      
      area_df_tmp = area_matrix[row_codes,]
      chems_df_tmp = cmp_matrix[row_codes,]
      probs_df_tmp = probs_matrix[row_codes,]
      mass_df_tmp = mass_matrix[row_codes,]
      RT_df_tmp = RT_matrix[row_codes,]
    } else {
      area_df_tmp = area_matrix[rtBYmass_long %in% row_codes,]
      chems_df_tmp = cmp_matrix[rtBYmass_long %in% row_codes,]
      probs_df_tmp = probs_matrix[rtBYmass_long %in% row_codes,]
      mass_df_tmp = mass_matrix[rtBYmass_long %in% row_codes,]
      RT_df_tmp = RT_matrix[rtBYmass_long %in% row_codes,]
    }
    
    if(nrow(area_df_tmp) == 0){next}
    # mass_df_tmp[is.na(mass_df_tmp)] = 0
    area_df_tmp[is.na(area_df_tmp)] = 0
    probs_df_tmp[is.na(probs_df_tmp)] = 0
    mass_df_tmp[is.na(mass_df_tmp)] = 0
    RT_df_tmp[is.na(RT_df_tmp)] = 0
    
    # area_df_tmp[mass_df_tmp == "58.041864811"] = 0
    # chems_df_tmp[mass_df_tmp == "58.041864811"] = "<NA>"
    # probs_df_tmp[mass_df_tmp == "58.041864811"] = 0
    # RT_df_tmp[mass_df_tmp == "58.041864811"] = 0
    # mass_df_tmp[mass_df_tmp == "58.041864811"] = 0
    
    area_df_tmp_apply2 = lapply(area_df_tmp, as.numeric)
    probs_df_tmp_apply2 = lapply(probs_df_tmp, as.numeric)
    mass_df_tmp_apply2 = lapply(mass_df_tmp, as.numeric)
    RT_df_tmp_apply2 = lapply(RT_df_tmp, as.numeric)
    
    probs_df_tmp_apply = lapply(seq_along(probs_df_tmp_apply2),function(i) {probs_df_tmp_apply2[[i]][RT_df_tmp_apply2[[i]] > 0]})
    probs_max_tmp_apply = lapply(seq_along(probs_df_tmp_apply2),function(i) {max(probs_df_tmp_apply2[[i]][RT_df_tmp_apply2[[i]] > 0])})
    
    mass_df_tmp_apply = lapply(seq_along(mass_df_tmp_apply2),function(i) {mass_df_tmp_apply2[[i]][RT_df_tmp_apply2[[i]] > 0]})
    RT_df_tmp_apply = lapply(seq_along(RT_df_tmp_apply2),function(i) {RT_df_tmp_apply2[[i]][RT_df_tmp_apply2[[i]] > 0]})
    
    best_RTs = lapply(seq_along(RT_df_tmp_apply),
                      function(i) {min(RT_df_tmp_apply[[i]][probs_df_tmp_apply[[i]] >= probs_max_tmp_apply[[i]] - 20])})
    previous_RTs = c(previous_RTs, best_RTs)
    best_RT = mean(as.vector(unlist(best_RTs)), na.rm = T)
    
    previous_RT_logical = as.matrix(unlist(lapply(seq_along(mass_df_tmp_apply2),function(i) {as.numeric(paste0(RT_df_tmp_apply2[[i]])) %in% previous_RTs})))
    if(any(previous_RT_logical)){
      area_df_tmp_apply = lapply(seq_along(area_df_tmp_apply2),function(i) {area_df_tmp_apply2[RT_df_tmp_apply2[[i]] %in% previous_RTs] = 0})
    }else{area_df_tmp_apply = area_df_tmp_apply2}     
    
    #RT_df_tmp_apply2 %in% previous_RTs)){area_df_tmp_apply2[RT_df_tmp_apply2 %in% previous_RTs] = 0}
    # area_df_tmp_apply2 = lapply(seq_along(area_df_tmp_apply2),function(i) {area_df_tmp_apply2[[i]][!area_df_tmp_apply2[[i]] %in% previous_areas]})
    # area_df_tmp_apply = lapply(seq_along(area_df_tmp_apply2),function(i) {area_df_tmp_apply2[[i]][RT_df_tmp_apply2[[i]] > 0 & between(RT_df_tmp_apply2[[i]], as.numeric(paste0(best_RT)), as.numeric(paste0(best_RT)) + RT_search_range)]})
    # previous_areas = c(previous_areas, as.numeric(paste0(area_df_tmp_apply)))
    # previous_RTs = c(previous_RTs, as.numeric(paste0(RT_df_tmp_apply)))
    
    
    best_masses = lapply(seq_along(mass_df_tmp_apply),
                         function(i) {mass_df_tmp_apply[[i]][probs_df_tmp_apply[[i]] == max(probs_df_tmp_apply[[i]], na.rm = T)]})
    
    
    # best_RTs = lapply(seq_along(RT_df_tmp_apply),
    #                   function(i) {RT_df_tmp_apply[[i]][probs_df_tmp_apply[[i]] >= max(probs_df_tmp_apply[[i]], na.rm = T)-15]})
    
    # best_low_RT = lapply(seq_along(best_RTs),function(i) {as.numeric(paste0(area_df_tmp_apply2[[i]][RT_df_tmp_apply2[[i]] == best_RTs[[i]]]))})
    
    # best_RT = as.numeric(paste0(best_low_RT)) #min(as.vector(unlist(best_RTs)), na.rm = T)
    all_best_RTs = c(all_best_RTs, best_RT)
    
    
    area_df_final_apply = lapply(seq_along(area_df_tmp_apply),function(i) {sum(as.numeric(paste0(area_df_tmp_apply2[[i]][RT_df_tmp_apply2[[i]] == best_RTs[[i]]])))})
    
    current_mass = as.numeric(paste0(masses_ordered[k]))
    if(current_mass == 999.999){
      best_RT = median(as.vector(unlist(best_RTs)), na.rm = T)
      current_mass = median(as.vector(unlist(best_masses)), na.rm = T)
    }
    # as.vector(unlist(area_df_tmp_apply))
    
    current_chem_row = cbind(as.numeric(paste0(best_RT)), 
                             as.numeric(paste0(current_mass)), 
                             as.character(paste0(chems_ordered[k])), 
                             rbind(as.vector(unlist(area_df_final_apply))))
    
    rownames(current_chem_row) = NULL
    colnames(current_chem_row) = c("RT",
                                   "Mass",
                                   "Chemical",
                                   colnames(area_matrix))

    final_list[[k]] = current_chem_row
  }
  final_df = do.call(rbind.data.frame,final_list)
  rownames(final_df) = NULL
  
  return(final_df)
}
