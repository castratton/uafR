#'@title exactoThese
#'
#'@description Takes categorated output as input and makes chemical subsets by the information in it.
#'
#'@details Provides a set of search chemicals for `mzExacto()`. User gets to select
#'chemicals based on information generated from `categorate()`.
#'
#'@param categoratedInput the lists output from `categorate()`
#'@param subsetBy specifies the list to subset by (Database, FMCS, or Library)
#'@param subsetArgs additional arguments to specify which values to subset by
#'@param subsetArgs2 additional arguments to specify which values to subset by
#'@param subset_input used when subsetting by FMCS information (e.g. molecular weight)
#'
#'@returns A vector of chemical names that meet each user-specification. Used as input
#'for `mzExacto()`.
#'
#'@examples
#'exactoThese(categoratedInput, subsetBy = "Database", subsetArgs = "All")
#'exactoThese(categoratedInput, subsetBy = "Database", subsetArgs = "reactives")
#'exactoThese(categoratedInput, subsetBy = "Database", subsetArgs = "LOTUS")
#'exactoThese(categoratedInput, subsetBy = "Database", subsetArgs = "KEGG")
#'exactoThese(categoratedInput, subsetBy = "Database", subsetArgs = "FEMA")
#'exactoThese(categoratedInput, subsetBy = "Database", subsetArgs = "FDA_SPL")
#'exactoThese(categoratedInput, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
#'exactoThese(categoratedInput, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Between", subset_input = c(125, 200))
#'@export

exactoThese = function(categoratedInput, subsetBy = "Database", subsetArgs = "All", subsetArgs2 = NA, subset_input = NA){
  subset_opts = c("Database", "FMCS", "Library")
  if(!(subsetBy %in% subset_opts)){stop("First Subset Argument Unrecognized!\n Your options are: 1. Database, 2. FMCS, 3. Library")}
  if(subsetBy == "Database"){
  dbSubset_args = list("All", "reactives", "LOTUS", "KEGG", "FEMA", "FDA_SPL",
                    c("reactives", "LOTUS"), c("reactives", "KEGG"),
                    c("reactives", "FEMA"), c("reactives", "FDA_SPL"),
                    c("LOTUS", "KEGG"), c("LOTUS", "FEMA"), c("LOTUS", "FDA_SPL"),
                    c("FEMA", "FDA_SPL"), c("reactives", "LOTUS", "KEGG"),
                    c("reactives", "LOTUS", "FEMA"), c("reactives", "LOTUS", "FDA_SPL"),
                    c("reactives", "KEGG", "FEMA"), c("reactives", "KEGG", "FDA_SPL"),
                    c("reactives", "FEMA", "FDA_SPL"), c("LOTUS", "KEGG", "FEMA"),
                    c("LOTUS", "KEGG", "FDA_SPL"), c("KEGG", "FEMA", "FDA_SPL"),
                    c("reactives", "LOTUS", "KEGG", "FEMA"),
                    c("reactives", "LOTUS", "KEGG", "FDA_SPL"),
                    c("LOTUS", "KEGG", "FEMA", "FDA_SPL"))
  dbSubset_expand = c("reactives", "LOTUS", "KEGG", "FEMA", "FDA_SPL")
  dbSubset_any_tmp1 = expand.grid(dbSubset_expand, dbSubset_expand)
  dbSubset_any_tmp2 = expand.grid(dbSubset_expand, dbSubset_expand, dbSubset_expand)
  dbSubset_any_tmp3 = expand.grid(dbSubset_expand, dbSubset_expand, dbSubset_expand, dbSubset_expand)
  dbSubset_any_tmp4 = expand.grid(dbSubset_expand, dbSubset_expand, dbSubset_expand, dbSubset_expand, dbSubset_expand)
  dbSubset_any = c(dbSubset_expand, do.call(paste, dbSubset_any_tmp1), do.call(paste, dbSubset_any_tmp2), do.call(paste, dbSubset_any_tmp3), do.call(paste, dbSubset_any_tmp4))

  reactives_logical = categoratedInput$Databases$reactives_df != "None"
  LOTUS_logical = categoratedInput$Databases$LOTUS_df != "None"
  KEGG_logical = categoratedInput$Databases$KEGG_df != "None"
  FEMA_logical = categoratedInput$Databases$FEMA_df != "None"
  FDA_SPL_logical = categoratedInput$Databases$FDA_SPL_df != "None"

  if(!(paste(subsetArgs, collapse = " ") %in% dbSubset_any)){stop("Database Subset Arguments Unrecognized, Please Try Again!\n If subsetting by multiple, their order matters: \n1. reactives, 2. LOTUS, 3. KEGG, 4. FEMA, 5. FDA_SPL")}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[1]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & LOTUS_logical & KEGG_logical & FEMA_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[2]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[3]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[LOTUS_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[4]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[KEGG_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[5]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[FEMA_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[6]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[7]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & LOTUS_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[8]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & KEGG_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[9]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & FEMA_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[10]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[11]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[LOTUS_logical & KEGG_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[12]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[LOTUS_logical & FEMA_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[13]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[LOTUS_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[14]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[FEMA_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[15]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & LOTUS_logical & KEGG_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[16]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & LOTUS_logical & FEMA_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[17]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & LOTUS_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[18]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & KEGG_logical & FEMA_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[19]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & KEGG_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[20]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & FEMA_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[21]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[LOTUS_logical & KEGG_logical & FEMA_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[22]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[LOTUS_logical & KEGG_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[23]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[KEGG_logical & FEMA_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[24]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & LOTUS_logical & KEGG_logical & FEMA_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[25]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[reactives_logical & LOTUS_logical & KEGG_logical & FDA_SPL_logical]}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[26]], collapse = " ")){exactoChems = categoratedInput$Databases$Chemical[LOTUS_logical & KEGG_logical & FEMA_logical & FDA_SPL_logical]}
 }

 if(subsetBy == "FMCS"){
  FMCS_subsetArgs = list("MW", "Rings", "Groups",
                         "Atoms", "NCharges")
  FMCS_subsetArgs2 = list("Equals", "Greater Than", "Less Than", "Between")
  FMCS_subset_val = subset_input


  if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) == FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) > FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) < FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$MW)) < FMCS_subset_val[2]]}

  if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) == FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) > FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) < FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$Rings)) < FMCS_subset_val[2]]}

  if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) == FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) > FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) < FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) < FMCS_subset_val[2]]}

  if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) == FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) > FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) < FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) < FMCS_subset_val[2]]}

  if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) == FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) > FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) < FMCS_subset_val]}
  if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$NCharges)) < FMCS_subset_val[2]]}
 }
 if(subsetBy == "Library"){
   library_groups = colnames(categoratedInput$FunctionalGroups)
   library_logical = library_groups %in% subsetArgs
   exactoChems = categoratedInput$FunctionalGroups$Chemical[categoratedInput$FunctionalGroups[,library_logical] != "No"]
 }
 tryCatch(return(unique(exactoChems[!is.na(exactoChems)])), error = function(error) {stop("Database Subset Arguments Unrecognized, Please Try Again!\n If subsetting by multiple, their order matters: \n1. reactives, 2. LOTUS, 3. KEGG, 4. FEMA, 5. FDA_SPL")})
}
