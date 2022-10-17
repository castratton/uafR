#'Post-merge function to quantify compound emission rates relative to internal or external standard(s).
#'
#'@description If using internal standard (IS), will use the 'standardify()' function with user-specified 
#'inputs. If using external standard (ES), will require an input "matrix" from which standard curves will 
#'be derived.
#'
#'@param input uafR 'theMerger()' output. 
#'

# data_in = unknowns_merged

standardifyIt = function(data_in, standard_type = "Internal", 
                         IS_ng = 190.5, IS_uL = 1, 
                         collect_time = 1, sample_amt = 1,
                         ES_calibration = NA){
 standardify = function(u, v, 
                        w, x, 
                        y, z){
  standardised = (u/v*w)/x/y/z
 }
 options(scipen = n)
 type_internal = standard_type == "Internal"
 type_external = standard_type == "External"
 
 size_sample_amt = length(sample_amt)
 
 data_clms = data_in$Area[,-c(1,2,3)]
 meta_clms = data_in$Area[,c(1,2,3)]
 IS_quants = as.numeric(paste0(data_in$IS[,-c(1,2,3)]))
 
 data_out_standard = data.frame(matrix(nrow = nrow(data_clms), ncol = 0))
 
 if (type_internal & size_sample_amt > 1){
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

   # for(i in 1:ncol(data_clms)){
   #   Input = data_clms[,i]
   #   IS = IS_quants[,i]
   #   sample_quant = sample_amt[,i]
   #   tmp = standardify(u = Input,
   #                     v = IS,
   #                     w = IS_ng,
   #                     x = IS_uL,
   #                     y = collect_time,
   #                     z = sample_quant)
   #   data_out_standard = data.frame(cbind(data_out_standard,
   #                                        tmp))
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
   
   # for(i in 1:ncol(data_clms)){
   #   Input = data_clms[,i]
   #   IS = IS_quants[,i]
   #   
   #   tmp = standardify(u = Input, 
   #                     v = IS, 
   #                     w = IS_ng, 
   #                     x = IS_uL, 
   #                     y = collect_time, 
   #                     z = sample_amt)
   #   data_out_standard = data.frame(cbind(data_out_standard, 
   #                                        tmp))
 } else {}
 data_out_standard = data.frame(cbind(meta_clms, data_out_standard))
 colnames(data_out_standard) = colnames(data_in$Area)
 return(data_out_standard)
 options(scipen = 0)
 }


 
#  if (type_internal){
#    
#    if (length(sample_amt) > 1){
#      for(i in 1:ncol(data_clms)){
#        Input = data_clms[,i]
#        IS = IS_quants[,i]
#        sample_quant = sample_amt[,i]
#        tmp = standardify(u = as.numeric(paste0(Input)), 
#                          v = as.numeric(paste0(IS)), 
#                          w = as.numeric(paste0(IS_ng)), 
#                          x = as.numeric(paste0(IS_uL)), 
#                          y = as.numeric(paste0(collect_time)), 
#                          z = as.numeric(paste0(sample_quant)))
#        data_out_standard = data.frame(cbind(data_out_standard, 
#                                      tmp))
#      }
#    } else {
#      for(i in 1:ncol(data_clms)){
#        Input = data_clms[,i]
#        IS = IS_quants[i]
#        
#        tmp = standardify(u = as.numeric(paste0(Input)), 
#                          v = as.numeric(paste0(IS)), 
#                          w = as.numeric(paste0(IS_ng)), 
#                          x = as.numeric(paste0(IS_uL)), 
#                          y = as.numeric(paste0(collect_time)), 
#                          z = as.numeric(paste0(sample_amt)))
#        data_out_standard = data.frame(cbind(data_out_standard, 
#                                             tmp))
#      }
#    }
#    
#  } else {}
#  
#  if (type_external){
#    
#    if (hasArg(ES_calibration)){
#      
#      ES_Abundance = ES_calibration$Component.Area
#      ES_ng = ES_calibration$Quantity
#      
#      model_list = list()
#      log_model = lm(ES_Abundance~log(ES_ng))
#      exponent_model = lm(ES_Abundance~exp(ES_ng))
#      linear_model = lm(ES_Abundance~ES_ng)
#      model_list[[length(model_list)+1]] = log_model
#      model_list[[length(model_list)+1]] = exponent_model
#      model_list[[length(model_list)+1]] = linear_model
#      
#      R_sq_set = c()
#      for(model in 1:length(model_list)){
#        model_current = model_list[[model]]
#        summary_current = summary(model_current)
#        R_sq_tmp = as.numeric(paste0(summary_current$adj.r.squared))
#        R_sq_set = c(R_sq_set, R_sq_tmp)
#      }
#      
#      best_fit = max(R_sq_set)
#      best_model = R_sq_set == best_fit
#      calibrating_curve = model_list[best_model]
#      eqn_coefficients = as.numeric(paste0(calibrating_curve[[1]]$coefficients))
#      
#      x_step = nrow(data_clms)
#      y_step = ncol(data_clms)
#      
#      area_standardized = matrix(nrow = x_step, ncol = y_step)
#      
#      for(x in 1:x_step){
#        for(y in 1:y_step){
#          current_area = as.numeric(paste0(data_clms[x,y]))
#          if(current_area == 0 | is.na(current_area)) next
#          area_standardized[x,y] = model_fun_use(current_area,
#                                                 eqn_coefficients[1],
#                                                 eqn_coefficients[2])
#        }
#        }
#      
#      if (length(sample_amt) > 1){
#        for(i in 1:ncol(data_clms)){
#          Input = data_clms[,i]
#          IS = IS_quants[,i]
#          sample_quant = sample_amt[,i]
#          tmp = standardify(u = Input, 
#                            v = IS, 
#                            w = IS_ng, 
#                            x = IS_uL, 
#                            y = collect_time, 
#                            z = sample_quant)
#          data_out_standard = data.frame(cbind(data_out_standard, 
#                                               tmp))
#        }
#        
#        
#      } else {
#        
#        for(i in 1:ncol(data_clms)){
#          Input = data_clms[,i]
#          IS = IS_quants[,i]
#          
#          tmp = standardify(u = Input, 
#                            v = IS, 
#                            w = IS_ng, 
#                            x = IS_uL, 
#                            y = collect_time, 
#                            z = sample_amt)
#          data_out_standard = data.frame(cbind(data_out_standard, 
#                                               tmp))
#        }
#        
#        
#    } else {print("Properly formatted calibration data not found!")}
#  } else {}
#  
#  if (!type_internal & !type_external){
#    print("Invalid entry!")} else {}
#  }
#  data_out_standard = data.frame(cbind(meta_clms, data_out_standard))
#  colnames(data_out_standard) = colnames(data_in$Area)
#  options(scipen = 0)
#  return(data_out_standard)
# }




 # if(standard_type == "Internal"){
 #    data_in_standard = data.frame(matrix(nrow = 0, ncol = ncol(data_in$Area)-3))
 #    for(i in 3:ncol(data_in$Area)){
 #       Input = data_in$Area[,i]
 #       IS = data_in$IS[,i]
 #       if(IS > 0){
 #          tmp = standardify(u = Input, 
 #                            v = IS, 
 #                            w = IS_ng, 
 #                            x = IS_uL, 
 #                            y = collect_time, 
 #                            z = sample_amt)
 #          data_in_standard = data.frame(data_in_standard, 
 #                                        tmp)
 #       } else {
          # IS = sum(data_in[,i-1][data_in$cmps_correct == int_standard])
          # tmp = standardify(u = Input, 
          #                   v = IS, 
          #                   w = IS_molecules, 
          #                   x = IS_ng, 
          #                   y = collect_time, 
          #                   z = sample_amt)
          # data_in_standard = data.frame(data_in_standard, 
          #                               tmp)
 #       }
 #    }
 # } else {} 
 # if(standard_type == "External" & hasArg(ES_calibration)){
 #       ES_Abundance = ES_calibration$Component.Area
 #       ES_ng = ES_calibration$Quantity
 #       
 #       # max_degree = length(unique(ES_ng))
 #       # model_list = list()
 #       # for(polynomial in 1:max_degree){
 #       #    model_tmp = lm(ES_Abundance~poly(ES_ng, polynomial, raw = T))
 #       #    model_list[[polynomial]] = model_tmp
 #       # }
 #       log_model = lm(ES_Abundance~log(ES_ng))
 #       # model_list[[length(model_list)+1]] = log_model
 #       
 #       # R_sq_set = c()
 #       # for(model in 1:length(model_list)){
 #       #    model_current = model_list[[model]]
 #       #    summary_current = summary(model_current)
 #       #    R_sq_tmp = as.numeric(paste0(summary_current$adj.r.squared))
 #       #    R_sq_set = c(R_sq_set, R_sq_tmp)
 #       # }
 #       # best_fit = max(R_sq_set)
 #       # best_model = R_sq_set == best_fit
 #       # calibrating_curve = model_list[best_model]
 #       # eqn_coefficients = as.numeric(paste0(calibrating_curve[[1]]$coefficients))
 #       
 #       # poly_fun1 = function(x, a1, b){
 #       #  y = a1*x + b
 #       #  return(y)
 #       # }
 #       # 
 #       # poly_fun2 = function(x, a1, a2, b){
 #       #  y = a2*x^2 + a1*x + b
 #       #  return(y)
 #       # }
 #       # 
 #       # poly_fun3 = function(x, a1, a2, a3, b){
 #       #  y = a3*x^3 + a2*x^2 + a1*x + b
 #       #  return(y)
 #       # }
 #       # 
 #       # poly_fun4 = function(x, a1, a2, a3, a4, b){
 #       #  y = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + b
 #       #  return(y)
 #       # }
 #       # 
 #       # poly_fun5 = function(x, a1, a2, a3, a4, a5, b){
 #       #  y = a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + b
 #       #  return(y)
 #       # }
 #       # 
 #       # poly_fun6 = function(x, a1, a2, a3, b){
 #       #  y = a6*x^6 + a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + b
 #       #  return(y)
 #       # }
 #       # 
 #       # poly_fun7 = function(x, a1, a2, a3, b){
 #       #  y = a7*x^7 + a6*x^6 + a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + b
 #       #  return(y)
 #       # }
 #       # 
 #       # poly_fun8 = function(x, a1, a2, a3, b){
 #       #  y = a8*x^8 + a7*x^7 + a6*x^6 + a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + b
 #       #  return(y)
 #       # }
 #       # poly_fun_list = list(poly_fun1, poly_fun2, poly_fun3, poly_fun4, poly_fun5, poly_fun6, poly_fun7, poly_fun8)
 #       # multiple possible functions depending on number of coefficients.
 #       poly_fun_use = poly_fun_list[[length(eqn_coefficients)-1]]
 #       area_to_STD = unknowns_merged$Area[,-c(1,2,3)]
 #       x_step = nrow(area_to_STD)
 #       y_step = ncol(area_to_STD)
 #       
 #       area_standardized = matrix(nrow = x_step, ncol = y_step)
 #       
 #       # if(length(eqn_coefficients) == 2){
 #       #    for(x in 1:x_step){
 #       #       for(y in 1:y_step){
 #       #          current_area = as.numeric(paste0(area_to_STD[x,y]))
 #       #          if(current_area == 0 | is.na(current_area)) next
 #       #          area_standardized[x,y] = poly_fun_use(current_area, 
 #       #                                                eqn_coefficients[1], 
 #       #                                                eqn_coefficients[2])
 #       #       }
 #       #    }
 #       # } else {}
 #       # 
 #       # if(length(eqn_coefficients) == 3){
 #       #    for(x in 1:x_step){
 #       #       for(y in 1:y_step){
 #       #          current_area = as.numeric(paste0(area_to_STD[x,y]))
 #       #          if(current_area == 0 | is.na(current_area)) next
 #       #          area_standardized[x,y] = poly_fun_use(current_area, 
 #       #                                                eqn_coefficients[1], 
 #       #                                                eqn_coefficients[2], 
 #       #                                                eqn_coefficients[3])
 #       #       }
 #       #    }
 #       # } else {}
 #       # 
 #       # if(length(eqn_coefficients) == 4){
 #       #    for(x in 1:x_step){
 #       #       for(y in 1:y_step){
 #       #          current_area = as.numeric(paste0(area_to_STD[x,y]))
 #       #          if(current_area == 0 | is.na(current_area)) next
 #       #          area_standardized[x,y] = poly_fun_use(current_area, 
 #       #                                                eqn_coefficients[1], 
 #       #                                                eqn_coefficients[2], 
 #       #                                                eqn_coefficients[3], 
 #       #                                                eqn_coefficients[4])
 #       #       }
 #       #    }
 #       # } else {}
 #       # 
 #       # if(length(eqn_coefficients) == 5){
 #       #    for(x in 1:x_step){
 #       #       for(y in 1:y_step){
 #       #          current_area = as.numeric(paste0(area_to_STD[x,y]))
 #       #          if(current_area == 0 | is.na(current_area)) next
 #       #          area_standardized[x,y] = poly_fun_use(current_area, 
 #       #                                                eqn_coefficients[1], 
 #       #                                                eqn_coefficients[2], 
 #       #                                                eqn_coefficients[3], 
 #       #                                                eqn_coefficients[4], 
 #       #                                                eqn_coefficients[5])
 #       #       }
 #       #    }
 #       # } else {}
 #       # 
 #       # if(length(eqn_coefficients) == 6){
 #       #    for(x in 1:x_step){
 #       #       for(y in 1:y_step){
 #       #          current_area = as.numeric(paste0(area_to_STD[x,y]))
 #       #          if(current_area == 0 | is.na(current_area)) next
 #       #          area_standardized[x,y] = poly_fun_use(current_area, 
 #       #                                                eqn_coefficients[1], 
 #       #                                                eqn_coefficients[2], 
 #       #                                                eqn_coefficients[3], 
 #       #                                                eqn_coefficients[4], 
 #       #                                                eqn_coefficients[5], 
 #       #                                                eqn_coefficients[6])
 #       #       }
 #       #    }
 #       # } else {}
 #       # 
 #       # if(length(eqn_coefficients) == 7){
 #       #    for(x in 1:x_step){
 #       #       for(y in 1:y_step){
 #       #          current_area = as.numeric(paste0(area_to_STD[x,y]))
 #       #          if(current_area == 0 | is.na(current_area)) next
 #       #          area_standardized[x,y] = poly_fun_use(current_area, 
 #       #                                                eqn_coefficients[1], 
 #       #                                                eqn_coefficients[2], 
 #       #                                                eqn_coefficients[3], 
 #       #                                                eqn_coefficients[4], 
 #       #                                                eqn_coefficients[5], 
 #       #                                                eqn_coefficients[6], 
 #       #                                                eqn_coefficients[7])
 #       #       }
 #       #    }
 #       # }else {}
 #       # 
 #       # if(length(eqn_coefficients) == 8){
 #       #    for(x in 1:x_step){
 #       #       for(y in 1:y_step){
 #       #          current_area = as.numeric(paste0(area_to_STD[x,y]))
 #       #          if(current_area == 0 | is.na(current_area)) next
 #       #          area_standardized[x,y] = poly_fun_use(current_area, 
 #       #                                                eqn_coefficients[1], 
 #       #                                                eqn_coefficients[2], 
 #       #                                                eqn_coefficients[3], 
 #       #                                                eqn_coefficients[4], 
 #       #                                                eqn_coefficients[5], 
 #       #                                                eqn_coefficients[6], 
 #       #                                                eqn_coefficients[7], 
 #       #                                                eqn_coefficients[8])
 #       #       }
 #       # } else {}
 #    } else {}
 # data_in$Area
  #  
  #  data_in_standard = area_standardized
  # }
 # data_in_standard = data.frame(cbind(, data_in_standard))
 # colnames(data_in_standard) = colnames(data_in$Area)
 # return(data_in_standard)
# }

