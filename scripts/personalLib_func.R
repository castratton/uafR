#'Formatting user's personal library
#'
#'@description
#' 'personalLib' takes a data set (formatted long or wide) and
#' returns named vectors of the provided data either as separate
#'  vectors or a list of vectors.
#'  
#'@details
#' This function allows easy manipulation of an imported data set 
#' for creation into "personal libraries" to use in downstream 
#' functions. Imported data set(s) can be in either long 
#' (containing one column for possible variable types and one 
#' column for the values of those variable types) or wide 
#' (containing a column labeled for each variable type with the 
#' values of that variable type in the corresponding rows) format. 
#' User may designate the format of the input data and the desired 
#' format for the output data. Data output options will either be 
#' in multiple separate character vectors (one per variable type) 
#' or a list (named "librarylist") containing the multiple separate
#' character vectors.
#' 
#' @param data A dataframe or matrix
#' @param input_Format either "long" or "wide"
#' @param output_Format either "vectors" or "list"
#' @returns multiple separate character vectors #' (one per 
#' variable type) when using ouput option "vectors"
#' @returns a list named "librarylist" containing the multiple 
#' separate character vectors when using output option "list" 
#' @examples personalLib(data, long, list)
#' @examples personalLib(data, wide, list)
#' @examples personalLib(data, long, vectors)
#' @examples personalLib(data, wide, vectors)
 
personalLib=function(data,input_Format,output_Format){
  if(input_Format == "long" && output_Format =="list"){
    df.2=na.omit(data)
    m=df.2[,c("compound","type")]
    dfwide=unstack(m)
    assign("librarylist",dfwide, envir = parent.frame())
    
  }
  
  if(input_Format == "wide" && output_Format =="list"){
    librarylist = list()
    for (i in 1:ncol(data)) {
      modify=data[[i]] 
      librarylist[[i]]=na.omit(modify)
    }
    name_lib=colnames(data)
    names(librarylist) <- name_lib
    assign("librarylist",librarylist, envir = parent.frame())
    
  }
  
  if(input_Format == "long" && output_Format =="vectors"){
    df.2=na.omit(data)
    m=df.2[,c("compound","type")]
    dfwide=unstack(m)
    for (i in 1:length(dfwide)){
      assign(names(dfwide[i]), dfwide[[i]], envir = parent.frame())
    }}
  
  if(input_Format == "wide" && output_Format =="vectors"){
    for (i in 1:ncol(data)){
      assign(colnames(data[i]), na.omit(data[[i]]), envir = parent.frame())
    }
  }
}
