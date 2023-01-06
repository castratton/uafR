#'Post-merge function to quantify compound emission rates relative to internal or external standard(s).
#'
#'@description If using internal standard (IS), will use the 'standardify()' function with user-specified 
#'inputs. If using external standard (ES), will require an input "matrix" from which standard curves will 
#'be derived.
#'
#'@param input uafR 'theMerger()' output. 
#'


ISmissing = function(data_in){
 scorify = function(u, v, w, x, y, z){
  standardised = (u/v*w)/x/y/z
 }
 
 data_clms = data_in$Area[,-c(1,2,3)]
 meta_clms = data_in$Area[,c(1,2,3)]
 
 data_out_standard = data.frame(matrix(nrow = nrow(data_clms), ncol = 0))
 
 for(i in 1:ncol(data_clms)){
  Input = data_clms[,i]
  IS = IS_quants[,i]
  sample_quant = sample_amt[,i]
  tmp = standardify(u = as.numeric(paste0(Input)), 
                    v = as.numeric(paste0(IS)), 
                    w = as.numeric(paste0(IS_ng)), 
                    x = as.numeric(paste0(IS_uL)), 
                    y = as.numeric(paste0(collect_time)), 
                    z = as.numeric(paste0(sample_quant)))
  data_out_standard = data.frame(cbind(data_out_standard, 
                                       tmp))
  }
 } else {}
 
 data_out_standard = data.frame(cbind(meta_clms, data_out_standard))
 colnames(data_out_standard) = colnames(data_in$Area)
 return(data_out_standard)
 options(scipen = 0)
}





