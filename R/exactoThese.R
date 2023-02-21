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
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "All")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "reactives")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "LOTUS")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "KEGG")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "FEMA")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "FDA_SPL")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
#'exactoThese(standard_categorated, subsetBy = "FMCS", subsetArgs = "MW",
#'subsetArgs2 = "Between", subset_input = c(125, 200))
#'@export

exactoThese = function(categoratedInput, subsetBy = "Database", subsetArgs = "All", subsetArgs2 = NA, subset_input = NA){
  subset_opts = c("Database", "FMCS", "Library")
  if(!(subsetBy %in% subset_opts)){stop("First Subset Argument Unrecognized!\n Your options are: 1. Database, 2. FMCS, 3. Library")}
  if(subsetBy == "Database"){
  dbSubset_args = list("All", "reactives", "LOTUS", "KEGG", "FEMA", "FDA_SPL",
                    c("reactives", "LOTUS"), c("reactives", "KEGG"),
                    c("reactives", "FEMA"), c("reactives", "FDA_SPL"),
                    c("LOTUS", "KEGG"), c("LOTUS", "FEMA"), c("LOTUS", "FDA_SPL"),
                    c("KEGG","FEMA"),c("KEGG","FDA_SPL"),
                    c("FEMA", "FDA_SPL"), c("reactives", "LOTUS", "KEGG"),
                    c("reactives", "LOTUS", "FEMA"), c("reactives", "LOTUS", "FDA_SPL"),
                    c("reactives", "KEGG", "FEMA"), c("reactives", "KEGG", "FDA_SPL"),
                    c("reactives", "FEMA", "FDA_SPL"), c("LOTUS", "KEGG", "FEMA"),
                    c("LOTUS", "KEGG", "FDA_SPL"), c("KEGG", "FEMA", "FDA_SPL"),
                    c("reactives", "LOTUS", "KEGG", "FEMA"),
                    c("reactives", "LOTUS", "KEGG", "FDA_SPL"),
                    c("LOTUS", "KEGG", "FEMA", "FDA_SPL"))

  dbSubset_expand = c("All", "reactives", "LOTUS", "KEGG", "FEMA", "FDA_SPL")
  dbSubset_any_tmp1 = expand.grid(dbSubset_expand, dbSubset_expand)
  dbSubset_any_tmp2 = expand.grid(dbSubset_expand, dbSubset_expand, dbSubset_expand)
  dbSubset_any_tmp3 = expand.grid(dbSubset_expand, dbSubset_expand, dbSubset_expand, dbSubset_expand)
  dbSubset_any_tmp4 = expand.grid(dbSubset_expand, dbSubset_expand, dbSubset_expand, dbSubset_expand, dbSubset_expand)
  dbSubset_any = c(dbSubset_expand, do.call(paste, dbSubset_any_tmp1), do.call(paste, dbSubset_any_tmp2), do.call(paste, dbSubset_any_tmp3), do.call(paste, dbSubset_any_tmp4))

  reactives_chems = categoratedInput$reactives$Chemical[categoratedInput$reactives$reactives != "None"]
  LOTUS_chems = categoratedInput$LOTUS$Chemical[categoratedInput$LOTUS$LOTUS != "None"]
  KEGG_chems = categoratedInput$KEGG$Chemical[categoratedInput$KEGG$KEGG != "None"]
  FEMA_chems = categoratedInput$FEMA$Chemical[categoratedInput$FEMA$FEMA != "None"]
  FDA_SPL_chems = categoratedInput$FDA_SPL$Chemical[categoratedInput$FDA_SPL$FDA_SPL != "None"]

  if(!(paste(subsetArgs, collapse = " ") %in% dbSubset_any)){stop("Database Subset Arguments Unrecognized, Please Try Again!\n If subsetting by multiple, their order matters: \n1. reactives, 2. LOTUS, 3. KEGG, 4. FEMA, 5. FDA_SPL")}
  if(paste(subsetArgs, collapse = " ") == "All"){exactoChems = Reduce(intersect, list(reactives_chems, LOTUS_chems))}
  if(paste(subsetArgs, collapse = " ") == "reactives"){exactoChems = unique(c(reactives_chems))}
  if(paste(subsetArgs, collapse = " ") == "LOTUS"){exactoChems = unique(c(LOTUS_chems))}
  if(paste(subsetArgs, collapse = " ") == "KEGG"){exactoChems = unique(c(KEGG_chems))}
  if(paste(subsetArgs, collapse = " ") == "FEMA"){exactoChems = unique(c(FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == "FDA_SPL"){exactoChems = unique(c(FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[7]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, LOTUS_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[8]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, KEGG_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[9]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[10]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[11]], collapse = " ")){exactoChems = Reduce(intersect, list(LOTUS_chems, KEGG_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[12]], collapse = " ")){exactoChems = Reduce(intersect, list(LOTUS_chems, FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[13]], collapse = " ")){exactoChems = Reduce(intersect, list(LOTUS_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[14]], collapse = " ")){exactoChems = Reduce(intersect, list(KEGG_chems,FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[15]], collapse = " ")){exactoChems = Reduce(intersect, list(KEGG_chems,FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[16]], collapse = " ")){exactoChems = Reduce(intersect, list(FEMA_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[17]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, LOTUS_chems, KEGG_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[18]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, LOTUS_chems, FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[19]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, LOTUS_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[20]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, KEGG_chems, FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[21]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, KEGG_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[22]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, FEMA_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[23]], collapse = " ")){exactoChems = Reduce(intersect, list(LOTUS_chems, KEGG_chems, FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[24]], collapse = " ")){exactoChems = Reduce(intersect, list(LOTUS_chems, KEGG_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[25]], collapse = " ")){exactoChems = Reduce(intersect, list(KEGG_chems, FEMA_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[26]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, LOTUS_chems, KEGG_chems, FEMA_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[27]], collapse = " ")){exactoChems = Reduce(intersect, list(reactives_chems, LOTUS_chems, KEGG_chems, FDA_SPL_chems))}
  if(paste(subsetArgs, collapse = " ") == paste(dbSubset_args[[28]], collapse = " ")){exactoChems = Reduce(intersect, list(LOTUS_chems, KEGG_chems, FEMA_chems, FDA_SPL_chems))}
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
