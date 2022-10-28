#'Prepares input GC/MS data for 'theMerger().'
#'
#'@description Every instance of every chemical across all samples is
#'placed on the same matrix with no overlaps.
#'
#'@param input uafR output, or csv containing relevant GC/MS data. 
#'

spreadOut = function(input){
 gcms_spread_area = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 gcms_spread_RT = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 gcms_spread_prob = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 gcms_spread_mz = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 gcms_spread_compound = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 gcms_spread_sample = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 gcms_spread_rtBYmass = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 gcms_spread_mass = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
 
 row.names(gcms_spread_area) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_area) = unique(input$File.Name)
 
 row.names(gcms_spread_RT) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_RT) = unique(input$File.Name)
 
 row.names(gcms_spread_prob) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_prob) = unique(input$File.Name)
 
 row.names(gcms_spread_mz) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_mz) = unique(input$File.Name)
 
 row.names(gcms_spread_compound) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_compound) = unique(input$File.Name)
 
 row.names(gcms_spread_sample) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_sample) = unique(input$File.Name)
 
 row.names(gcms_spread_rtBYmass) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_rtBYmass) = unique(input$File.Name)
 
 row.names(gcms_spread_mass) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
 colnames(gcms_spread_mass) = unique(input$File.Name)
 ##################################################################################################
 tentative_identities = as.character(paste0(input$Compound.Name))
 num_tentative = length(tentative_identities)
 tentative_masses = c()
 
 for (chem in 1:num_tentative){
   current_CMP = tentative_identities[chem]
   chem_cid = get_cid(tentative_identities[chem])
   chem_cid = paste0(chem_cid[[1,2]])
   # chem_cid = gsub("NA", 
   #                    "180", 
   #               chem_cid)
   # chem_cid = gsub("[c\\\"() ]",
   #                    "",
   #               chem_cid)
   # chem_cid = gsub("\n", 
   #                    "", 
   #               chem_cid)
   # chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
   #                 chem_cid,
   #                   '/SDF')
   
   # chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid,'/JSON?heading=GC-MS')
   
   # chem_info = tryCatch(jsonlite::fromJSON(chem_url), 
   #                      error = function(error) {return("None")})
   if(chem_cid == "NA"){chem_cid = 180}
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
   
   step_printer = c(1:num_tentative)
   step_cmper = chem/25
   
   if(step_cmper %in% step_printer | chem == 1){
     cat(paste0('\n', '[Current/Total]', ' |--Exact Mass--|--CMP Name 1--|--CMP Name 2--|', '\n'))
   } else {}
   
   names(mass_tmp) = NULL
   names(CMP_name_tmp) = NULL
   all_tmp = paste0(as.numeric(paste0(mass_tmp)), 
                    ';', 
                    CMP_name_tmp, 
                    ';', 
                    current_CMP)
   all_print = paste0(mass_tmp, 
                      ' | ', 
                      CMP_name_tmp, 
                      ' | ',
                      current_CMP)
   cat(paste0('[', chem, ' / ', num_tentative, ']', ': ', all_print, '\n'))
   tentative_masses = c(tentative_masses, as.numeric(paste0(mass_tmp)))
 }
 
 
 # tentative_mass_split = strsplit(tentative_masses, ';')
 # tentative_mass = as.numeric(paste0(sapply(tentative_mass_split, '[[', 1)))
 # tentative_CMP1 = as.character(paste0(sapply(tentative_mass_split, '[[', 2)))
 # tentative_CMP2 = as.character(paste0(sapply(tentative_mass_split, '[[', 3)))
 # tentative_df = data.frame(cbind(tentative_mass, tentative_CMP1, tentative_CMP2))
 # colnames(tentative_df) = c("Mass", "Chemical1", "Chemical2")
 # rownames(tentative_df) = NULL
 
 # tentative_df[sort(tentative_df$Mass, decreasing = T),]
 # tentative_mass_df = do.call(rbind.data.frame, tentative_mass_split)
 # colnames(tentative_mass_df) = NULL
 ################################################################################
 
 sample_names = unique(input$File.Name)
 compound_names = unique(input$Compound.Name)
 compound_area = unique(input$Component.Area)
 compound_RT = unique(input$Component.RT)
 compound_match = unique(input$Match.Factor)
 print("Spreading out your data, one datum at a time.")
 
 for (j in 1:nrow(gcms_spread_area)){
  area_tmp = paste0(input$Component.Area[j])
  compound_tmp = paste0(input$Compound.Name[j])
  sample_tmp = paste0(input$File.Name[j])
  mz_tmp = paste0(input$Base.Peak.MZ[j])
  prob_tmp = paste0(input$Match.Factor[j])
  RT_tmp = paste0(input$Component.RT[j])
  mass_tmp = paste0(tentative_masses[j])
  rtBYmass_tmp = paste0(RT_tmp,' | ', mass_tmp)
  
  gcms_spread_area[rownames(gcms_spread_area) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_area) == paste0(sample_tmp)] = area_tmp
  gcms_spread_compound[rownames(gcms_spread_compound) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_compound) == paste0(sample_tmp)] = compound_tmp
  gcms_spread_sample[rownames(gcms_spread_sample) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_sample) == paste0(sample_tmp)] = sample_tmp
  gcms_spread_mz[rownames(gcms_spread_mz) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_mz) == paste0(sample_tmp)] = mz_tmp
  gcms_spread_prob[rownames(gcms_spread_prob) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_prob) == paste0(sample_tmp)] = prob_tmp
  gcms_spread_RT[rownames(gcms_spread_RT) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_RT) == paste0(sample_tmp)] = RT_tmp
  gcms_spread_rtBYmass[rownames(gcms_spread_rtBYmass) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_rtBYmass) == paste0(sample_tmp)] = rtBYmass_tmp
  gcms_spread_mass[rownames(gcms_spread_mass) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_mass) == paste0(sample_tmp)] = mass_tmp
 }
 
 print("Losing rownames (data ID codes).")
 row.names(gcms_spread_area) = NULL
 row.names(gcms_spread_RT) = NULL
 row.names(gcms_spread_prob) = NULL
 row.names(gcms_spread_mz) = NULL
 row.names(gcms_spread_compound) = NULL
 row.names(gcms_spread_sample) = NULL
 row.names(gcms_spread_mass) = NULL
 row.names(gcms_spread_rtBYmass) = NULL
 
 gcms_list = list(gcms_spread_area, gcms_spread_compound, 
                  gcms_spread_mz, gcms_spread_prob, 
                  gcms_spread_RT, gcms_spread_mass, 
                  gcms_spread_rtBYmass)
 names(gcms_list) = c("Area", "Compounds", 
                      "MZ", "MatchFactor", 
                      "RT", "Mass", 
                      "rtBYmass")
 return(gcms_list)
}
