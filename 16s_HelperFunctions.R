#16S_HelperFunctions.R

# Source required functions and load required libraries.

source("R_helper_functions.R")

# Helper functions.

check_phyloseqObject <- function(phyloseqObject, taxRank = "Phylum") {
  # Performs sanity checks on a phyloseq object.
  #
  # Args:
  #   phyloseqObject: An object of class phyloseq.
  #   taxRank: The taxonomic rank at which to retrieve unique taxa.
  #
  # Returns:
  #   The sample variables, number of taxa, rank names, and unique taxa of the
  #     phyloseq object.
  
  # Error handling 
  # Check that the phyloseq library is loaded
  if (!require("phyloseq", character.only = TRUE)) {
    stop("Please make sure the \"phyloseq\" library is loaded.")
  }
  # Check that a phyloseq object has been passed in
  if ((class(phyloseqObject)) != "phyloseq") {
    stop("Please check that you called this function on an object of type 
         \"phyloseq.\"")
  }
  
  # Get the rank names in the phyloseq object. Check that a user-specified rank
  #   name is available in the phyloseq object.
  rankNames <- rank_names(phyloseqObject)
  if ((taxRank %in% rankNames) == "FALSE") {
    warning("Warning: the specified taxonomic rank was not one of the 
            following: ", 
                  paste(as.character(rankNames), collapse = ", "))
    warning(". Setting taxonomic rank to \"Phylum.\"")
    taxRank = "Phylum"
  }
  
  # Perform sanity checks and write out to user
  sampleVars <- sample_variables(phyloseqObject)
  taxaNumber <- ntaxa(phyloseqObject)
  uniqueTaxa <- get_taxa_unique(phyloseqObject, taxRank)
  
  writeLines(paste("\nSample Variables:\n\n", paste(as.character(sampleVars), 
                                                    collapse = ", ")))
  writeLines(paste("\nTaxa Number:", as.character(taxaNumber)))
  writeLines(paste("\nRank Names:", paste(as.character(rankNames), 
                                          collapse = ", ")))
  writeLines(paste("\nChecked for unique taxa at the ", 
                   as.character(taxRank), "level."))
  writeLines(paste("\nUnique Taxa: ", paste(as.character(uniqueTaxa), 
                                            collapse = ", ")))
  
}