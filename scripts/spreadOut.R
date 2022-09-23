#'Prepares input GC/MS data for 'theMerger().'
#'
#'@description Every instance of every chemical across all samples is
#'placed on the same matrix with no overlaps.
#'
#'@param input uafR output, or csv containing relevant GC/MS data. 
#'

spreadOut = function(input){
 gcms_spread_area = data.frame(matrix(ncol = length(unique(input$Sample.Name)), nrow = length(input$Component.RT)))
 gcms_spread_RT = data.frame(matrix(ncol = length(unique(input$Sample.Name)), nrow = length(input$Component.RT)))
 gcms_spread_prob = data.frame(matrix(ncol = length(unique(input$Sample.Name)), nrow = length(input$Component.RT)))
 gcms_spread_mz = data.frame(matrix(ncol = length(unique(input$Sample.Name)), nrow = length(input$Component.RT)))
 gcms_spread_compound = data.frame(matrix(ncol = length(unique(input$Sample.Name)), nrow = length(input$Component.RT)))
 gcms_spread_sample = data.frame(matrix(ncol = length(unique(input$Sample.Name)), nrow = length(input$Component.RT)))
 
 row.names(gcms_spread_area) = paste0(input$Component.RT, input$Base.Peak.MZ)
 colnames(gcms_spread_area) = unique(input$Sample.Name)
 
 row.names(gcms_spread_RT) = paste0(input$Component.RT, input$Base.Peak.MZ)
 colnames(gcms_spread_RT) = unique(input$Sample.Name)
 
 row.names(gcms_spread_prob) = paste0(input$Component.RT, input$Base.Peak.MZ)
 colnames(gcms_spread_prob) = unique(input$Sample.Name)
 
 row.names(gcms_spread_mz) = paste0(input$Component.RT, input$Base.Peak.MZ)
 colnames(gcms_spread_mz) = unique(input$Sample.Name)
 
 row.names(gcms_spread_compound) = paste0(input$Component.RT, input$Base.Peak.MZ)
 colnames(gcms_spread_compound) = unique(input$Sample.Name)
 
 row.names(gcms_spread_sample) = paste0(input$Component.RT, input$Base.Peak.MZ)
 colnames(gcms_spread_sample) = unique(input$Sample.Name)
 
 sample_names = unique(input$Sample.Name)
 compound_names = unique(input$Compound.Name)
 compound_area = unique(input$Component.Area)
 compound_RT = unique(input$Component.RT)
 compound_match = unique(input$Match.Factor)
 print("Spreading out your data, one datum at a time.")
 for (j in 1:nrow(gcms_spread_area)){
  area_tmp = paste0(input$Component.Area[j])
  compound_tmp = paste0(input$Compound.Name[j])
  sample_tmp = paste0(input$Sample.Name[j])
  mz_tmp = paste0(input$Base.Peak.MZ[j])
  prob_tmp = paste0(input$Match.Factor[j])
  RT_tmp = paste0(input$Component.RT[j])
  
  gcms_spread_area[rownames(gcms_spread_area) == paste0(RT_tmp, mz_tmp), colnames(gcms_spread_area) == paste0(sample_tmp)] = area_tmp
  gcms_spread_compound[rownames(gcms_spread_compound) == paste0(RT_tmp, mz_tmp), colnames(gcms_spread_compound) == paste0(sample_tmp)] = compound_tmp
  gcms_spread_sample[rownames(gcms_spread_sample) == paste0(RT_tmp, mz_tmp), colnames(gcms_spread_sample) == paste0(sample_tmp)] = sample_tmp
  gcms_spread_mz[rownames(gcms_spread_mz) == paste0(RT_tmp, mz_tmp), colnames(gcms_spread_mz) == paste0(sample_tmp)] = mz_tmp
  gcms_spread_prob[rownames(gcms_spread_prob) == paste0(RT_tmp, mz_tmp), colnames(gcms_spread_prob) == paste0(sample_tmp)] = prob_tmp
  gcms_spread_RT[rownames(gcms_spread_RT) == paste0(RT_tmp, mz_tmp), colnames(gcms_spread_RT) == paste0(sample_tmp)] = RT_tmp
 }
 print("Lose rownames (data ID codes).")
 row.names(gcms_spread_area) = NULL
 row.names(gcms_spread_RT) = NULL
 row.names(gcms_spread_prob) = NULL
 row.names(gcms_spread_mz) = NULL
 row.names(gcms_spread_compound) = NULL
 row.names(gcms_spread_sample) = NULL
 
 gcms_list = list(gcms_spread_area, gcms_spread_compound, 
                  gcms_spread_mz, gcms_spread_prob, 
                  gcms_spread_RT)
 names(gcms_list) = c("Area", "Compounds", 
                      "MZ", "MatchFactor", 
                      "RT")
 return(gcms_list)
}
