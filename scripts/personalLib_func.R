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
