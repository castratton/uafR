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

mzMerger <- function(data_in){
 area_matrix = data_in$Area  
 probs_matrix = data_in$MatchFactor
 cmp_matrix = data_in$Compounds
 RT_matrix = data_in$RT
 mz_matrix = data_in$MZ
 
 CMPs_long = c()
 
 for (z in 1:nrow(cmp_matrix)){
  CMPs_tmp = cmp_matrix[z,][!is.na(cmp_matrix[z,])]
  CMPs_long = c(CMPs_long, CMPs_tmp)
 }
 
 chemicals = c(unique(as.character(CMPs_long)))
 
 # chemical_list = list()
 chemical_mz_all = c()
 
 area_decays = data_in$Area
 probs_decays = data_in$MatchFactor
 cmp_decays = data_in$Compounds
 RT_decays = data_in$RT
 mz_decays = data_in$MZ
 
 i = 4
 
 for (i in 1:length(chemicals)){
  chemical_mz_this = c()
  
  chemical_cid = get_cid(chemicals[i])
  
  if(!is.na(chemical_cid[,2])){
   current_chemical_cid = paste0(as.matrix(chemical_cid[,2]))
   current_chemical_name = paste0(as.matrix(chemical_cid[,1]))
  } else {next}
  
  # current_chemical = chemical_cid
  chemical_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', current_chemical_cid,'/JSON?heading=GC-MS')
  chemical_info = tryCatch(jsonlite::fromJSON(chemical_url), 
                              error = function(error) {return("None")})
  
  if(chemical_info == "None"){next}
  
  chemical_json_dat = data.frame(unlist(chemical_info$Record$Section))
  colnames(chemical_json_dat) = "info"
  mz_rows = grep("[[:digit:]]+\\.[[:digit:]]\\ [[:digit:]]+\\.[[:digit:]]", chemical_json_dat$info)
   
  chemical_mass_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
                       current_chemical_cid,
                       '/SDF')
  chemical_mass_info = tryCatch(read.SDFset(chemical_mass_url), 
                          error = function(error) {return("None")})
  chemical_mass_tmp = chemical_mass_info@SDF[[1]]@datablock[names(chemical_mass_info@SDF[[1]]@datablock) == "PUBCHEM_EXACT_MASS"]
  chemical_name_tmp = chemical_mass_info@SDF[[1]]@datablock[names(chemical_mass_info@SDF[[1]]@datablock) == "PUBCHEM_IUPAC_CAS_NAME"]
  
  if(length(mz_rows) == 0){
    backup_search = c("Section.Section.Information.Value.Number2",
                      "Section.Section.Information.Value.Number3",
                      "Section.Section.Information.Value.Number4",
                      "Section.Section.Information.Value.Number5")
    
    mz_rows_backup = paste0(chemical_json_dat$info[row.names(chemical_json_dat) %in% backup_search])
    
    mz_rows = mz_rows_backup
    if(length(mz_rows) == 0){next}
    chemical_mz_matches = unique(paste0(unlist(strsplit(paste0(chemical_json_dat[mz_rows,])," "))))
    chemical_mz_matches2 = gsub("\\.0\\>","",chemical_mz_matches)
    
    backup_switch = T
  } else {
    chemical_mz_matches = unique(paste0(unlist(strsplit(paste0(chemical_json_dat[mz_rows,])," "))))
    chemical_mz_matches2 = gsub("\\.0\\>","",chemical_mz_matches)
  }
  # chemical_list = list(chemical_list, chemical_mz_matches2)
  # names(chemical_list)[i] = chemical_names[i]
  chemical_mz_this = c(chemical_mz_this, chemical_mz_matches2)
  
  mz_count = 1:length(chemical_mz_this)
  odds = mz_count[mz_count%%2 != 0]
  evens = mz_count[mz_count%%2 == 0]
  
  mz_primary = chemical_mz_this[odds]
  mz_secondary = chemical_mz_this[evens]
  
  chemical_mz_this_odd = mz_primary
  # chemical_mz_all = mz_primary
  
  mz_print1 = paste(mz_primary, collapse = " | ")
  mz_print2 = paste(mz_secondary, collapse = " | ")
  # chemical_mz_all=NULL
  
  if(backup_switch){
    chemical_mz_this = chemical_mz_this[!chemical_mz_this %in% chemical_mz_all]
  } else {
    chemical_mz_this = chemical_mz_this_odd[!chemical_mz_this_odd %in% chemical_mz_all]
  }
  
  if(length(chemical_mz_this)<1){next}
  print(paste0("[", i, "/", length(chemicals), "] ", "Merging ",paste0(chemical_name_tmp)," by top m/z peaks: ", paste(chemical_mz_this, collapse = " | ")))
  
  mz_rows = c()
  for (mz in 1:length(chemical_mz_this)){
   mz_rows_tmp = which(mz_matrix == chemical_mz_this[mz], arr.ind = T)
   mz_rows_tmp = paste0(mz_rows_tmp[,1])
   mz_rows = c(mz_rows, mz_rows_tmp)
  }
  mz_rows = unique(mz_rows)
  
  chemical_mz_all = c(chemical_mz_all, mz_primary)
  
  area_matrix_chem = area_decays[mz_rows,]
  probs_matrix_chem = probs_decays[mz_rows,]
  cmp_matrix_chem = cmp_decays[mz_rows,]
  RT_matrix_chem = RT_decays[mz_rows,]
  mz_matrix_chem = mz_decays[mz_rows,]
  
  cmp_row = c()
  probs_row = c()
  RT_row = c()
  mz_row = c()
  
  chemical_row_list = list()
  # rownames(chemical_row_list) = NULL
  # colnames(chemical_row_list) = c("RT", 
  #                            "Mass", 
  #                            "Chemical", 
  #                            colnames(data_in$Area))
  
  area_matrix_chem[is.na(area_matrix_chem)] = 0
  probs_matrix_chem[is.na(probs_matrix_chem)] = 0
  # cmp_matrix_chem[is.na(cmp_matrix_chem)] = 0
  RT_matrix_chem[is.na(RT_matrix_chem)] = 0
  mz_matrix_chem[is.na(mz_matrix_chem)] = 0
  
  if (nrow(area_matrix_chem) > 0){
  
   for(c in 1:ncol(RT_matrix_chem)){
    
    # cmp_row_tmp = as.character(paste0(cmp_matrix_chem[,c]))
    probs_row_tmp = as.numeric(paste0(probs_matrix_chem[,c]))
    max_prob_1 = max(probs_row_tmp) 
    
    RT_row_tmp = as.numeric(paste0(RT_matrix_chem[,c]))
    RT_row_tmp_1 = RT_row_tmp[probs_row_tmp == max_prob_1] 
    mz_row_tmp = as.numeric(paste0(mz_matrix_chem[,c]))
    mz_row_tmp_1 = mz_row_tmp[probs_row_tmp == max_prob_1]
    
    probs_row = c(probs_row, max_prob_1)
    RT_row = c(RT_row, RT_row_tmp_1)
    mz_row = c(mz_row, mz_row_tmp_1)
   }
   
   max_prob = max(probs_row)
   
   best_identified = length(probs_row[probs_row > max_prob - 10])
   best_RT = mean(RT_row[probs_row %in% tail(sort(probs_row), best_identified)])
   # best_area = mean(STD_area_1[STD_probs %in% tail(sort(STD_probs),best_identified)])
   # mz_Standard[RT_Standard <= best_RT - 3 | RT_Standard >= best_RT + 3]
   area_matrix_chem_tmp = area_matrix_chem[RT_row >= best_RT - 3.5 | RT_row <= best_RT + 3.5,]
   # area_Standard = area_Standard[area_Standard >= best_area/10 | area_Standard <= best_area*10,]
   area_matrix_chem = do.call(cbind.data.frame, area_matrix_chem_tmp)
   area_matrix_chem[is.na(area_matrix_chem)] = 0
   
   probs_matrix_chem_tmp = probs_matrix_chem[RT_row >= best_RT - 3.5 | RT_row <= best_RT + 3.5,]
   # probs_Standard = probs_Standard[probs_Standard >= best_probs/10 | probs_Standard <= best_probs*10,]
   probs_matrix_chem = do.call(cbind.data.frame, probs_matrix_chem_tmp)
   probs_matrix_chem[is.na(probs_matrix_chem)] = 0
   
   area_row = c()
   
   
   for (col in 1:ncol(area_matrix_chem)){
    area_row_tmp = as.numeric(paste0(area_matrix_chem[,col]))
    summed_area = sum(area_row_tmp[area_row_tmp > 0])
    
    area_row = c(area_row, summed_area)
   }
   
   chemical_row = data.frame(cbind(best_RT, chemical_mass_tmp, chemical_name_tmp, rbind(area_row)))
   
   rownames(chemical_row) = NULL
   colnames(chemical_row) = c("RT", 
                         "Mass", 
                         "Chemical", 
                         colnames(data_in$Area))
   chemical_row_list[[length(chemical_row_list) + 1]] = chemical_row
   # chemical_df[i,] = chemical_row
  } else {next}
  
  area_decays = area_decays[!rownames(area_decays) %in% mz_rows,]
  probs_decays = probs_decays[!rownames(probs_decays) %in% mz_rows,]
  cmp_decays = cmp_decays[!rownames(cmp_decays) %in% mz_rows,]
  RT_decays = RT_decays[!rownames(RT_decays) %in% mz_rows,]
  mz_decays = mz_decays[!rownames(mz_decays) %in% mz_rows,]
  
  row.names(area_decays) = NULL
  row.names(probs_decays) = NULL
  row.names(cmp_decays) = NULL
  row.names(RT_decays) = NULL
  row.names(mz_decays) = NULL
  
 }
 chemical_df = do.call(rbind.data.frame, chemical_row_list)
 row.names(chemical_df) = NULL
 colnames(chemical_df) = c("RT", 
                           "Mass", 
                           "Chemical", 
                           colnames(data_in$Area))
 
 return(chemical_df)
}
