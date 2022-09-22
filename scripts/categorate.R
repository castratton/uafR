#'Categorizing chemical input based on published and proprietary 
#'information
#'
#'@description Uses 'get_cid()' to convert chemical names into unique 
#'compound identifier values and pull records containing info from 5 
#'databases accessed through PubChem. Unmatched compounds are removed. 
#'It also uses 'dictionate()' to categorize chemicals by common, 
#'functionally-defined structures. Output is a list of data frames with 
#'the output for the $databases, $FMCS, and $FunGroups data generation.
#'
#'@param input A vector containing chemical names in IUPAC notation 
#'(preferred).
#'

categorate = function(compounds){
  compound_cids = get_cid(compounds)
  compound_cids = compound_cids[!is.na(compound_cids$cid),]
  
  compound_SDFs = suppressWarnings(dictionate(make.unique(compounds)))
  print("3-D chemical structure data acquired from PubChem.")
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
  
  print("Scouring the internet for information on these chemicals.")
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
    if (length(reactives_final) < 1){
      reactives_final = "None"
    } else {}
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
    if (length(LOTUS_final) < 1){
      LOTUS_final = "None"
    } else {}
    LOTUS_df = LOTUS_final
    
    if (KEGG_info != "None"){
      KEGG_value = as.matrix(unlist(KEGG_info))
    } else {KEGG_value = "None"}
    if (!is.na(KEGG_value[2])){
      KEGG_final = KEGG_value[substr(KEGG_value,1,5)=="KEGG:"]
    } else {KEGG_final = "None"}
    if (length(KEGG_final) < 1){
      KEGG_final = "None"
    } else {}
    KEGG_df = KEGG_final
    
    if (FEMA_info != "None"){
      FEMA_value = as.matrix(unlist(FEMA_info))
    } else {FEMA_value = "None"}
    if (!is.na(FEMA_value[2])){
      FEMA_final = FEMA_value[row.names(FEMA_value)=="Record.Section.Section.Information.Value.StringWithMarkup.String"]
    } else {FEMA_final = "None"}
    if (length(FEMA_final) < 1){
      FEMA_final = "None"
    } else {}
    FEMA_df = unlist(strsplit(FEMA_final,","))
    
    if (MeSH_info != "None"){
      MeSH_value = as.matrix(unlist(MeSH_info))
    } else {MeSH_value = "None"}
    if (!is.na(MeSH_value[2])){
      MeSH_final = MeSH_value[row.names(MeSH_value)=="Record.Section.Section.Information.Value.StringWithMarkup.String"]
    } else {MeSH_final = "None"}
    if (length(MeSH_final) < 1){
      MeSH_final = "None"
    } else {}
    MeSH_df = MeSH_final
    
    if (FDA_SPL_info != "None"){
      FDA_SPL_value = as.matrix(unlist(FDA_SPL_info))
    } else {FDA_SPL_value = "None"}
    if (!is.na(FDA_SPL_value[2])){
      FDA_SPL_final = FDA_SPL_value[row.names(FDA_SPL_value)=="Record.Reference.Name"]
    } else {FDA_SPL_final = "None"}
    if (length(FDA_SPL_final) < 1){
      FDA_SPL_final = "None"
    } else {}
    FDA_SPL_df = FDA_SPL_final
    
    CMP_info_row = as.data.frame(t(rbind.data.frame(reactives_df, LOTUS_df, KEGG_df, FEMA_df, FDA_SPL_df, Chemical)))
    colnames(CMP_info_row) = Chem_data_source
    CMP_info_df = rbind(CMP_info_df, CMP_info_row)
    row.names(CMP_info_df) = NULL
  }
  print("Done with the internet, moving on to chemical structures.")
  
  SDF_info_df = data.frame(matrix(ncol = 9, nrow = 0))
  SDF_columns = c("Chemical", "MW", "MF", "Rings", 
                  "Groups", "GroupCounts", "Atom", 
                  "AtomCounts", "NCharges")
  
  colnames(SDF_info_df) = SDF_columns
  cid(compound_SDFs) <- makeUnique(cid(compound_SDFs))
  # SDF = 322
  for(SDF in 1:length(compound_SDFs)){
    
    Ncharges = tryCatch(sapply(bonds(compound_SDFs[SDF], type="charge"),length), error = function(error) {return("None")})
    atom_counts = tryCatch(atomcountMA(compound_SDFs[SDF], addH = FALSE), error = function(error) {return("None")})
    group_counts = tryCatch(groups(compound_SDFs[SDF], type = "countMA"), error = function(error) {return("None")})
    ring_counts = tryCatch(rings(compound_SDFs[SDF], upper = 6, type = "count", arom = TRUE, type = "countMA"), error = function(error) {return("None")})
    Chemical = tryCatch(compound_SDFs[[SDF]][[4]][[9]], error = function(error) {return("None")})
    MF = tryCatch(MF(compound_SDFs[SDF], addH=FALSE), error = function(error) {return("None")})
    MW = tryCatch(MW(compound_SDFs[SDF], addH=FALSE), error = function(error) {return("None")})
    
    atom_count_atoms = tryCatch(colnames(atom_counts), error = function(error) {return("None")})
    group_count_groups = tryCatch(rownames(as.matrix(unlist(group_counts))), error = function(error) {return("None")})
    if(is.null(group_count_groups)){
      group_count_groups = "None"
    }
    print("Chemicals dissected, will now match 3D structures.")
    SDF_info_row = as.data.frame(t(rbind.data.frame(as.vector(Chemical), 
                                                    as.vector(MW), as.vector(MF), 
                                                    as.vector(ring_counts),
                                                    group_count_groups,
                                                    as.vector(group_counts), 
                                                    atom_count_atoms,
                                                    as.vector(atom_counts), 
                                                    as.vector(Ncharges))))
    batch_test_set = tryCatch(fmcsBatch(compound_SDFs[SDF], functional_SDFs, au = 0, bu = 0), error = function(error) {return(batch_test_set = rbind(rep(0, 5), rep(0, 5)))})
    
    if (max(as.vector(batch_test_set[,5])) > 0.95)
    {
      functional_match = names(batch_test_set[,5][batch_test_set[,5] > 0.95])
      match_SDFs = functional_SDFs[functional_SDFs@ID %in% functional_match]
      functional_match = as.integer(gsub("CMP", "", functional_match))
      functional_row = data.frame(t(rep("No", ncol(functional_df))))
      colnames(functional_row) = functional_identities
      functional_row$Chemical = Chemical
      # batch_test_set[,5][names(batch_test_set[,5]) %in% functional_match]
      functional_row[names(batch_test_set[,5]) %in% functional_match] = "Yes"
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
