#'Accessing chemical structure files from PubChem in SDF format
#'
#'@description Uses 'get_cid()' to convert chemical names into unique 
#'compound identifier values and pull records from PubChem. Unmatched 
#'compounds are replaced with a filler compound (acetone) to maintain 
#'the original length of the input. Automatically bypasses PubChem 
#'limits on URL length and maximum download time from a server. 
#'
#'@param input A vector containing a list of chemical names.
#'

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
