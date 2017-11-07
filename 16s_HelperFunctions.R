#16S_HelperFunctions.R

###########################################################
# Source required functions and load required libraries. #
###########################################################

source("R_helper_functions.R")
#load_R_libraries("ArgumentCheck")

#####################
# Helper Functions #
#####################

check_phyloseqObject <- function(phyloseqObject, taxRank = "Phylum") {

  # Perform Checks: 
  
  # Check that the phyloseq library is loaded into the calling script's environment.  Don't try to proceed otherwise.
  if (!require("phyloseq", character.only = TRUE)) {
    stop("Please make sure the \"phyloseq\" library is loaded.")
  }
  
  # Check that a phyloseq object has been passed in.  Don't try to proceed otherwise.
  if ((class(phyloseqObject)) != "phyloseq") {
    stop("Please check that you called this function on an object of type \"phyloseq.\"")
  }
  
  # Get the rank names in the phyloseq object and use that to check that the taxRank argument specified is ok.
  rankNames <- rank_names(phyloseqObject)
  if ((taxRank %in% rankNames) == "FALSE") {
    warning("Warning: the specified taxonomic rank was not one of the following: ", 
                  paste(as.character(rankNames), collapse = ", "))
    warning(". Setting taxonomic rank to \"Phylum.\"")
    taxRank = "Phylum"
  }
  
  # Check the phyloseq object.
  sampleVars <- sample_variables(phyloseqObject)
  taxaNumber <- ntaxa(phyloseqObject)
  uniqueTaxa <- get_taxa_unique(phyloseqObject, taxRank)
  
  writeLines(paste("\nSample Variables:\n\n", paste(as.character(sampleVars), collapse = ", ")))
  writeLines(paste("\nTaxa Number:", as.character(taxaNumber)))
  writeLines(paste("\nRank Names:", paste(as.character(rankNames), collapse = ", ")))
  writeLines(paste("\nChecked for unique taxa at the ", as.character(taxRank), "level."))
  writeLines(paste("\nUnique Taxa: ", paste(as.character(uniqueTaxa), collapse = ", ")))
  
}