#'Post Alignment Toolbox
#'
#'@description Function that allows the user of the package to take
#'control of their data by subsetting according to user defined 
#'parameters. 
#'
#'Decontaminate, a function to eliminate compounds that are not of
#'interest from your dataset based on a user-defined character string option 
#'for exact or non-exact matches, default is exact.

# data_in = unknowns_spread

decontaminate <- function(data_in){
 contaminants = c("Perfluorotributylamine","Methanol","Benzene",
                  "Toluene","M-Xylene",
                  "O-Xylene","P-Xylene",
                  "Trichloroethane","dimethoxy(dimethyl)silane")
 
 area_matrix = data_in$Area  
 probs_matrix = data_in$MatchFactor
 cmp_matrix = data_in$Compounds
 RT_matrix = data_in$RT
 mz_matrix = data_in$MZ
 
 contaminant_cids = get_cid(contaminants)
 contaminant_cids = paste0(as.matrix(contaminant_cids[,2]))
 
 contaminant_mz_all = c()
 
 for (i in 1:length(contaminant_cids)){
  current_contaminant = contaminant_cids[i]
  contaminant_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', current_contaminant,'/JSON?heading=GC-MS')
  
  contaminant_info = tryCatch(jsonlite::fromJSON(contaminant_url), 
                     error = function(error) {return("None")})
  contaminant_json_dat = data.frame(unlist(contaminant_info$Record$Section))
  colnames(contaminant_json_dat) = "info"
  mz_rows = grep("[[:digit:]]+\\.[[:digit:]]\\ [[:digit:]]+\\.[[:digit:]]", contaminant_json_dat$info)
  
  contaminant_mz_matches = unique(paste0(unlist(strsplit(paste0(contaminant_json_dat[mz_rows,])," "))))
  contaminant_mz_matches2 = gsub("\\.0\\>","",contaminant_mz_matches)
  contaminant_mz_all = c(contaminant_mz_all, contaminant_mz_matches2)
 }
 mz_rows = c()
 for (mz in 1:length(contaminant_mz_all)){
  mz_rows_tmp = which(mz_matrix == contaminant_mz_all[mz], arr.ind = T)
  mz_rows_tmp = paste0(mz_rows_tmp[,1])
  mz_rows = c(mz_rows, mz_rows_tmp)
  }
 mz_rows = unique(mz_rows)
 
 area_matrix = data_in$Area[!rownames(data_in$Area) %in% mz_rows,]
 probs_matrix = data_in$MatchFactor[!rownames(data_in$MatchFactor) %in% mz_rows,]
 cmp_matrix = data_in$Compounds[!rownames(data_in$Compounds) %in% mz_rows,]
 RT_matrix = data_in$RT[!rownames(data_in$RT) %in% mz_rows,]
 mz_matrix = data_in$MZ[!rownames(data_in$MZ) %in% mz_rows,]
 
 gcms_list = list(area_matrix, cmp_matrix, 
                  mz_matrix, probs_matrix, 
                  RT_matrix)
 names(gcms_list) = c("Area", "Compounds", 
                      "MZ", "MatchFactor", 
                      "RT")
 return(gcms_list)
}
