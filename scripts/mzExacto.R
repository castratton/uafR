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
# IS = "Tetradecane"
chemicals = c("Acetic ester", "Ethyl hexanoate", "Octanal", "Undecane", "Methyl salicylate")

mzExacto <- function(data_in, chemicals){
  area_matrix = data_in$Area
  probs_matrix = data_in$MatchFactor
  cmp_matrix = data_in$Compounds
  RT_matrix = data_in$RT
  mz_matrix = data_in$MZ
  mass_matrix = data_in$Mass
  rtBYmass_matrix = data_in$rtBYmass
  
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
    current_CMP = chemicals[chem]
    chem_cid = get_cid(chemicals[chem])
    chem_cid = paste0(chem_cid[[1,2]])
    
    chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid,'/JSON?heading=GC-MS')
    
    chem_info = tryCatch(jsonlite::fromJSON(chem_url), 
                         error = function(error) {return("None")})
    
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
    
    if (length(mz_rows_1) > 0) {
      mz_rows = mz_rows_1
    } else {
      if (length(mz_rows_2) > 0) {
        mz_rows = mz_rows_2
      } else {
        if(length(mz_rows_3) > 0) {
          mz_rows = mz_rows_3
        } else {next}
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
    chem_mz_matches2 = mz_primary[order(mz_secondary,decreasing = T)]
    
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
  names(mz_mass_list) = chemicals
  
  mass_each_chem = sapply(mz_mass_list, '[[', 1)
  mz_each_chem = sapply(mz_mass_list, '[[', -1)
  
  mass_orderer = as.numeric(paste0(mass_each_chem))

  masses_ordered = mass_each_chem[order(mass_orderer)]
  mz_mass_ordered = mz_each_chem[order(mass_orderer)]
  chems_ordered = names(masses_ordered)
  
  CMPs_long = c()
  RTs_long = c()
  masses_long = c()
  MZs_long = c()
  rtBYmass_long = c()
  MZs_published = list()
  
  for (z in 1:nrow(cmp_matrix)){
    CMPs_tmp = cmp_matrix[z,][!is.na(cmp_matrix[z,])]
    CMPs_long = c(CMPs_long, CMPs_tmp)
    
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
    
    chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
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
    
    if (length(mz_rows_1) > 0) {
      mz_rows = mz_rows_1
    } else {
      if (length(mz_rows_2) > 0) {
        mz_rows = mz_rows_2
      } else {
        if(length(mz_rows_3) > 0) {
          mz_rows = mz_rows_3
        } else {next}
      }
    }
    
    chem_mz_matches = unique(paste0(unlist(strsplit(paste0(chem_json_dat[mz_rows,])," "))))
    chem_mz_matches2 = gsub("\\.0\\>","",chem_mz_matches)
    chem_mz_matches2 = gsub("\\.[[:digit:]]0\\>","",chem_mz_matches2)
    
    mz_count = 1:length(chem_mz_matches2)
    odds = mz_count[mz_count%%2 != 0]
    evens = mz_count[mz_count%%2 == 0]
    
    mz_primary = chem_mz_matches2[odds]
    if(chem_cid == 180){mz_primary = NA}
    mz_secondary = chem_mz_matches2[evens]
    chem_mz_matches2 = mz_primary[order(mz_secondary, decreasing = T)]
    
    mz_print = paste(chem_mz_matches2, collapse = " | ")
    
    step_printer = c(1:num_unique_CMPs)
    step_cmper = z/7
    
    if(step_cmper %in% step_printer | z == 1){
      cat(paste0('\n', '[Current/Total]', ' |--Chemical Name--|--Top MZ Peaks--|', '\n'))
    } else {}
    
    names(mass_tmp) = NULL
    
    all_print = paste0(CMPs_tmp,
                       ': ',
                       mz_print)
    cat(paste0('[', z, '/', nrow(cmp_matrix), ']', '-', all_print, '\n'))
    
    MZs_published[[z]] = mz_primary
  }
  if(length(MZs_published) < length(masses_long)){MZs_published[[length(masses_long)]] = NA}
  names(MZs_published) = paste0(RTs_long, ' | ', masses_long)
  
  row_find_df = data.frame(cbind(CMPs_long, as.numeric(paste0(masses_long)), MZs_long))
  rownames(row_find_df) = NULL
  colnames(row_find_df) = c("Chemical", "Mass", "MZ")
  
  probs_found = c()
  RT_found = c()
  area_found = c()
  
  masses_ordered = mass_each_chem[order(mass_orderer)]
  mz_mass_ordered = mz_each_chem[order(mass_orderer)]
  chems_ordered = names(masses_ordered)
  k = 1
  for(k in 1:length(chems_ordered)){
    row_codes = c()
    if(k == 1){
      all_masses = as.numeric(paste0(masses_long))
      current_mass = as.numeric(paste0(masses_ordered[k]))
      
      MZ_possibles = MZs_published[all_masses <= current_mass]
    } else {
      all_masses = as.numeric(paste0(masses_long))
      previous_mass = as.numeric(paste0(masses_ordered[k-1]))
      current_mass = as.numeric(paste0(masses_ordered[k]))
      
      MZ_possibles = MZs_published[all_masses > previous_mass & all_masses <= current_mass]
    }
    # p = 3
    for(p in 1:length(MZ_possibles)){
      tryCatch(length(MZ_possibles[[p]]), error = function(error) {next})
      if(length(MZ_possibles[[p]]) == 0){next}
      bad_mass_list = strsplit(names(MZ_possibles[p]),' | ')
      bad_mass = as.numeric(paste0(bad_mass_list[[1]][3]))
      # bad_RT = as.numeric(paste0(bad_mass_list[[1]][1]))
      if(bad_mass == 58.04186){next}
      # if(bad_RT < current_RT - 3 | bad_RT > current_RT + 3){next}
      if(any(mz_mass_ordered[[k]] %in% MZ_possibles[[p]])){
        row_codes_tmp = names(MZ_possibles[MZ_possibles[[p]] %in% mz_mass_ordered[[k]]])
        row_codes = c(row_codes, row_codes_tmp)
      } else {next}
    }
    
    area_df_tmp = area_matrix[rtBYmass_long %in% row_codes,]
    chems_df_tmp = cmp_matrix[rtBYmass_long %in% row_codes,]
    probs_df_tmp = probs_matrix[rtBYmass_long %in% row_codes,]
    mass_df_tmp = mass_matrix[rtBYmass_long %in% row_codes,]
    RT_df_tmp = RT_matrix[rtBYmass_long %in% row_codes,]

    # area_found = c()
    best_probs = c()
    # area_df_tmp[is.na(area_df_tmp)] = 0
    probs_df_tmp[is.na(probs_df_tmp)] = 0
    RT_df_tmp[is.na(RT_df_tmp)] = 0
    for(c in 1:ncol(RT_df_tmp)){
      # area_tmp = as.numeric(paste0(area_df_tmp[,c]))
      # area_tmp = sum(area_tmp[area_tmp > 0])
      probs_tmp = as.numeric(paste0(probs_df_tmp[,c]))
      probs_tmp = max(probs_tmp)
      best_probs = c(best_probs, probs_tmp)
      RT_tmp = as.numeric(paste0(RT_df_tmp[,c]))
      RT_tmp = mean(RT_tmp[RT_tmp > 0])
      
      # area_found = c(area_found, area_tmp)
      # probs_found = c(probs_found, probs_tmp)
      RT_found = c(RT_found, RT_tmp)
    }
    
    # best_identified = max(probs_found)
    RT_new = c()
    RT_found2 = RT_df_tmp
    probs_found2 = probs_df_tmp
    max_prob = as.numeric(paste0(max(probs_found2[probs_found2>0])))
    min_prob = as.numeric(paste0(min(probs_found2[probs_found2>0])))
    probs_diff = max_prob - min_prob
    # probs_found2[!(probs_found2 %in% best_probs),] = 0
    # i = 1
    for(i in 1:length(best_probs)){
      RT_new_tmp = RT_found2[probs_found2 >= best_probs[i]-probs_diff]
      RT_new_tmp = as.numeric(paste0(RT_new_tmp))
      RT_new = c(RT_new, RT_new_tmp)
    }
    # logical_RT = unique(RT_found2)
    
    # RT_found2 = RT_found2[probs_found2[probs_found2 == best_probs]  & mass_df_tmp != 58.04186]
    # 
    # RT_found2 = as.numeric(paste0(RT_found2))
    
    # if(length(logical_RT) <= 1){
    #   RT_foundB = RT_df_tmp
    #   RT_foundB = RT_foundB[RT_foundB > 0 & mass_df_tmp != 58.04186]
    #   # mass_foundB = as.numeric(paste0(mass_df_tmp[!is.na(mass_df_tmp)]))
    #   
    #   RT_found2 = as.numeric(paste0(RT_foundB))
    # } else {
    #   RT_foundA = RT_df_tmp[chems_df_tmp == chems_ordered[k]]
    #   RT_foundA[is.na(RT_foundA)] = 0
    #   RT_foundA = RT_foundA[RT_foundA > 0 & mass_df_tmp != 58.04186]
    #   RT_found2 = as.numeric(paste0(RT_foundA))
    # }
    max_RT = as.numeric(paste0(max(RT_new[RT_new>0])))
    min_RT = as.numeric(paste0(min(RT_new[RT_new>0])))
    RT_diff = max_RT - min_RT
    best_RT = median(RT_new[!is.na(RT_new)])
    
    area_found = c()
    
    # as.data.frame(RT_df_tmp)
    
    RT_df_tmp[is.na(RT_df_tmp)] = 0
    best_RT_all = c()
    # RT_df_tmp2 = do.call(cbind.data.frame, RT_df_tmp)
    # RT_df_tmp2 = data.frame(RT_df_tmp2)
    for(col in 1:ncol(RT_df_tmp)){
      RT_tentative = as.numeric(paste0(RT_df_tmp[,col]))
      best_RT_tmp = RT_df_tmp[RT_tentative >= min_RT & RT_tentative <= best_RT, col]
      best_RT_all = c(best_RT_all, best_RT_tmp)
    }
    
    best_RT = mean(as.numeric(paste0(best_RT_all[best_RT_all != "0"])))
    
    RT_array_tmp = RT_matrix
    rownames(RT_array_tmp) = NULL
    RT_array_tmp[is.na(RT_array_tmp)] = 0
    # RT_set = rownames(RT_array_tmp[RT_array_tmp <= best_RT])
    RT_array_tmp[]
    
    area_df_tmp[is.na(area_df_tmp)] = 0
    RT_array = RT_df_tmp >= min_RT & RT_df_tmp <= best_RT
    
    # for(col_2 in 1:ncol(RT_array)){}
    best_RT_rows = as.numeric(paste0(names(RT_array[,1])))
    area_df_tmpB = area_matrix
    rownames(area_df_tmpB) = NULL
    area_df_tmp2 = area_df_tmpB[rownames(area_df_tmpB) %in% best_RT_rows,]
    area_df_tmp2[is.na(area_df_tmp2)] = 0
    
    probs_df_tmpB = probs_matrix
    rownames(probs_df_tmpB) = NULL
    probs_df_tmp2 = probs_df_tmpB[rownames(probs_df_tmpB) %in% best_RT_rows,]
    probs_df_tmp2[is.na(probs_df_tmp2)] = 0
    # probs_df_tmp2 = probs_df_tmp[best_RT_rows,]
    # probs_df_tmp2[is.na(probs_df_tmp2)] = 0
    
    RT_df_tmpB = RT_matrix
    rownames(RT_df_tmpB) = NULL
    RT_df_tmp2 = RT_df_tmpB[rownames(RT_df_tmpB) %in% best_RT_rows,]
    RT_df_tmp2[is.na(RT_df_tmp2)] = 0
    # RT_df_tmp2 = RT_df_tmp[best_RT_rows,]
    # RT_df_tmp2[is.na(RT_df_tmp2)] = 0
    
    for(c in 1:ncol(RT_df_tmp)){
      area_tmp2 = as.numeric(paste0(area_df_tmp2[,c]))
      area_tmp2 = sum(area_tmp2[area_tmp2 > 0])
      probs_tmp2 = as.numeric(paste0(probs_df_tmp2[,c]))
      probs_tmp2 = max(probs_tmp2)
      RT_tmp2 = as.numeric(paste0(RT_df_tmp2[,c]))
      RT_tmp2 = mean(RT_tmp2[RT_tmp2 > 0])
      
      area_found = c(area_found, area_tmp2)
      probs_found = c(probs_found, probs_tmp2)
      RT_found = c(RT_found, RT_tmp2)
    }
    
    current_chem_row = cbind(as.numeric(paste0(best_RT)), 
                             as.numeric(paste0(masses_ordered[k])), 
                             as.character(paste0(chems_ordered[k])), 
                             rbind(area_found))
    
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
