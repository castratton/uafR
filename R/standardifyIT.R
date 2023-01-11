#'@title standardifyIt
#'
#'@description Post-merge function to quantify compound emission rates relative to internal or external
#'standard(s). If using internal standard (IS), will use the contained `standardify()` function with
#'user-specified inputs. If using external standard (ES), an input "matrix" from which standard curves
#'can be derived.
#'
#'@param data_in mzExacto output
#'@param standard_type specifies type of standardization to perform ("Internal" or "External")
#'@param standard_used specifies the standard used as an internal standard (name your chemical)
#'@param IS_ng specifies the quantified number of molecules in the standard
#'@param IS_uL specifies the amount of standard added to the samples
#'@param collect_time specifies how long samples were collected for (e.g. hour(s))
#'@param sample_amt specifies how many individuals/items the samples were collected from
#'@param ES_calibration input matrix from which calibrations will be fit for external standardization
#'
#'@importFrom methods hasArg
#'@importFrom stats lm
#'@export

standardifyIt = function(data_in, standard_type = "Internal",
                         standard_used = "Tetradecane",
                         IS_ng = 190.5, IS_uL = 1,
                         collect_time = 1, sample_amt = 1,
                         ES_calibration = NA){
 standardify = function(u, v, w, x, y, z){standardised = (u/v*w)/x/y/z}
 model_fun_use = NULL
 options(scipen = 999)
 type_internal = standard_type == "Internal"
 type_external = standard_type == "External"

 size_sample_amt = length(sample_amt)

 data_clms = data_in[,-c(1,2,3,4)]
 meta_clms = data_in[,c(1,2,3,4)]

 if(type_internal){
    data_out_standard = data.frame(matrix(nrow = nrow(data_clms)-1, ncol = 0))
    IS = standard_used
    IS_quants = as.numeric(paste0(data_clms[meta_clms$Compound == IS,]))
    data_clms = data_clms[meta_clms$Compound != IS,]
    meta_clms = meta_clms[meta_clms$Compound != IS,]
 }else{
    data_out_standard = data.frame(matrix(nrow = nrow(data_clms), ncol = 0))
 }



 if (type_internal & size_sample_amt > 1){
   for(i in 1:ncol(data_clms)){
      Input = data_clms[, i]
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

 if (type_internal & size_sample_amt == 1){
   for(i in 1:ncol(data_clms)){
      Input = data_clms[,i]
      IS = IS_quants[i]

      tmp = standardify(u = as.numeric(paste0(Input)),
                        v = as.numeric(paste0(IS)),
                        w = as.numeric(paste0(IS_ng)),
                        x = as.numeric(paste0(IS_uL)),
                        y = as.numeric(paste0(collect_time)),
                        z = as.numeric(paste0(sample_amt)))
      data_out_standard = data.frame(cbind(data_out_standard,
                                           tmp))
   }
 } else {}

 if (type_external & size_sample_amt > 1 & hasArg(ES_calibration)){
   ES_Abundance = ES_calibration$Component.Area
   ES_ng = ES_calibration$Quantity

   model_list = list()
   log_model = lm(ES_Abundance~log(ES_ng))
   exponent_model = lm(ES_Abundance~exp(ES_ng))
   linear_model = lm(ES_Abundance~ES_ng)
   model_list[[length(model_list)+1]] = log_model
   model_list[[length(model_list)+1]] = exponent_model
   model_list[[length(model_list)+1]] = linear_model

   R_sq_set = c()

   for(model in 1:length(model_list)){
     model_current = model_list[[model]]
     summary_current = summary(model_current)
     R_sq_tmp = as.numeric(paste0(summary_current$adj.r.squared))
     R_sq_set = c(R_sq_set, R_sq_tmp)
   }

   best_fit = max(R_sq_set)
   best_model = R_sq_set == best_fit
   calibrating_curve = model_list[best_model]
   eqn_coefficients = as.numeric(paste0(calibrating_curve[[1]]$coefficients))

   x_step = nrow(data_clms)
   y_step = ncol(data_clms)

   area_standardized = matrix(nrow = x_step, ncol = y_step)

   for(x in 1:x_step){
     for(y in 1:y_step){
       current_area = as.numeric(paste0(data_clms[x,y]))
       if(current_area == 0 | is.na(current_area)) next
       area_standardized[x,y] = model_fun_use(current_area,
                                              eqn_coefficients[1],
                                              eqn_coefficients[2])
     }
   }

 } else {}

 if (type_external & size_sample_amt == 1 & hasArg(ES_calibration)){
   ES_Abundance = ES_calibration$Component.Area
   ES_ng = ES_calibration$Quantity

   model_list = list()
   log_model = lm(ES_Abundance~log(ES_ng))
   exponent_model = lm(ES_Abundance~exp(ES_ng))
   linear_model = lm(ES_Abundance~ES_ng)
   model_list[[length(model_list)+1]] = log_model
   model_list[[length(model_list)+1]] = exponent_model
   model_list[[length(model_list)+1]] = linear_model

   R_sq_set = c()
   for(model in 1:length(model_list)){
     model_current = model_list[[model]]
     summary_current = summary(model_current)
     R_sq_tmp = as.numeric(paste0(summary_current$adj.r.squared))
     R_sq_set = c(R_sq_set, R_sq_tmp)
   }

   best_fit = max(R_sq_set)
   best_model = R_sq_set == best_fit
   calibrating_curve = model_list[best_model]
   eqn_coefficients = as.numeric(paste0(calibrating_curve[[1]]$coefficients))

   x_step = nrow(data_clms)
   y_step = ncol(data_clms)

   area_standardized = matrix(nrow = x_step, ncol = y_step)

   for(x in 1:x_step){
     for(y in 1:y_step){
       current_area = as.numeric(paste0(data_clms[x,y]))
       if(current_area == 0 | is.na(current_area)) next
       area_standardized[x,y] = model_fun_use(current_area,
                                              eqn_coefficients[1],
                                              eqn_coefficients[2])
     }
   }
 } else {}
 data_out_standard = data.frame(cbind(meta_clms, data_out_standard))
 colnames(data_out_standard) = colnames(data_in)
 return(data_out_standard)
 options(scipen = 0)
}