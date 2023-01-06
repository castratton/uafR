#'General compound matching (3D structure) function.
#'
#'@description Runs 'fmcsBatch()' on each individual compound in a 
#'set of compounds in SDF format (e.g. output of 'dictionate()'). 
#'Compares the input to a user-defined batch of chemicals also in 
#'SDF format. 
#'
#'@param SDF_input An SDF object containing the compounds that will 
#'be individually matched to a user-defined set. 
#'
#'@param comparison_set An SDF object that the input set will be 
#'matched against.
#'
#'@param target_value A number between 0 and 1 that specifies the 
#'portion of the molecule (in 3 dimensions) that must be identical 
#'to the comparison molecule for a match to be made. The calculated 
#'value is called the Tanimoto Coefficient.

TanimotoR = function(SDF_input, comparison_set, target_value){
 output = c()
 for(j in 1:length(SDF_input)){
  batch_test_set = fmcsBatch(SDF_input[[j]],comparison_set)
  # print(batch_test_set)
  if (max(as.vector(batch_test_set[,4])) > target_value)
  {
   output[j] = "Match"
  }else
  {
   output[j] = "No Match"
  }
 }
 return(output)
}