#'Import function for uafR

#'@description Imports data from the file location specified with 'path' and 
#'transforms it to the format required for GCAligR's align_chromatogram
#'algorithm. '.' and '-' in sample names and column names are replaced with 
#''_' to avoid issues with align_chromatograms.
#'The default format is a .csv file with the peak list of all samples 
#'that is output of Agilent's Masshunter Unknowns Analysis software. 
#'Alternatively, format = "folder" can be specified to compile the dataset from
#'individual peak list .csv files for each sample.
#'
#'@param path Character string specifying the file location of the data to be
#'imported
#'
#'@param format Character string with two options,".csv" indicating one .csv 
#'file containing all the input data, "folder" indicating that sample
#'chromatograms are in individual files in the folder that 'path' specifies.
#'Default is ".csv".
#'
#'@param ID_colname Name of the column containing the sample name. Relevant for
#'the ".csv" format parameter input. Default is "Sample_Name".
#'
#'@param subset_peaklist Specifies if subsetting according to the subsequent
#'parameters is desired. 
#'
#'@param rt_colname Name of column with retention time information.
#'
#'@param rtlow_cutoff Lower threshold for retention times to be included in
#'dataset.
#'@param rthigh_cutoff Upper threshold for retention times to be included in
#'dataset.
#'
#'@param matchf_colname Column name for the compound identification confidence 
#'value (Match factor).
#'
#'@param min_matchf Minimum required ID confidence for row to be retained.
#'
#'@param area_colname Column name for the area values of each peak.
#'
#'@param min_area Minimum area value for rows to be retained.
#'



gcms_data_in <- function(path, format = ".csv", ID_colname = "Sample.Name",
                         subset_peaklist= F,
                         rt_colname = "Component.RT", rtlow_cutoff = 0, 
                         rthigh_cutoff= 100, matchf_colname = "Match.Factor", 
                         min_matchf = 0, area_colname = "Component.Area",
                         min_area = 0){
  #here we could implent a checkpoint that checks if all the required variables
  #loaded and otherwise reports an error message
if (format == ".csv") {
  unknowns_all = read.csv(path)
  if (subset_peaklist == T) {
    unknowns_all <- subset(unknowns_all, unknowns_all[rt_colname] <= rthigh_cutoff & 
                             unknowns_all[rt_colname] >= rtlow_cutoff & 
                             unknowns_all[matchf_colname] >= min_matchf & 
                             unknowns_all[area_colname] >= min_area)
  }
  #Replacing '.' and '-' with '_' to avoid GCAlignR error message. Keep?
  names(unknowns_all)<- gsub('\\.','_',names(unknowns_all))
  names(unknowns_all)<- gsub('-','_',names(unknowns_all))
  ID_colname <- gsub('\\.','_',ID_colname)
  ID_colname <- gsub('-','_',ID_colname)
  #Splits dataframe into list format
  out <- split(unknowns_all,unknowns_all[[ID_colname]])
  names(out) <- gsub('\\.','_',names(out))
  names(out) <- gsub('-','_',names(out))
 }

#Alternative input from folder of .csv files, file names will be sample names!!
if (format == "folder") {
  fileList = list.files(path="C:/Intercropping/Data_Kernza", pattern=".csv", full.names = TRUE)
  GC_input = lapply(fileList, fread)
  #subset function, notyet tested for this type of import, please check before publishing
  if (subset_peaklist == T ) {
    for (i in 1:length(GC_input)) {
      GC_input[[i]] <- subset(GC_input[[i]], GC_input[[i]][rt_colname] <= rthigh_cutoff & 
                                GC_input[[i]][rt_colname] >= rtlow_cutoff & 
                                GC_input[[i]][matchf_colname] >= min_matchf & 
                                GC_input[[i]][area_colname] >= min_area)
    }
  }
  # GC_input is without list item labels, the next step imports labels from the file names
  fileNames = list.files(path="C:/Intercropping/Data_Kernza", pattern = NULL, all.files = TRUE, 
                         full.names = FALSE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE
  )
  #removes .csv
  fileNames = substr(fileNames,1,nchar(fileNames)-4)
  fileNames = data.frame(fileNames)
  fileNames = gsub('\\.','_', fileNames)
  fileNames = gsub('-','_', fileNames)
}
return(out)
}

#The following GCAlignR function could be integrated??
#check_input(out))