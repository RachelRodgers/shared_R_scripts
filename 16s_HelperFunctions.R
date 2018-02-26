#16S_HelperFunctions.R

################################################################################
# Colorblind-friendly pallets for plotting
cbPaletteGrey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPaletteBlack <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# For fills use
#   scale_fill_manual(values = cbPalletteGrey)

# For line and point colors use
#   scale_colour_manual(values = cbPalleteGrey)



################################################################################
# Helper functions.

CheckPhyloseqObject <- function(phyloseqObject, taxRank = "Phylum") {
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

RemoveMissingTaxa <- function(physeq) {
  # Removes any taxa that are missing (have a taxa sum of 0) from a physeq object
  # Error handling
  # Check that a phyloseq object has been passed in
  if ((class(physeq)) != "phyloseq") {
    stop("Please check that you called this function on an object of type 
         \"phyloseq.\"")
  }
  
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  #return(physeq)
}


MakeAbundanceDF <- function(physeq, tax.rank, abundance.filter = 0.01) {
  # Creates an abundance data frame at a given taxonomic rank for a phyloseq object
  # for abundance bar plots
  abundance.df <- physeq %>%
    tax_glom(taxrank = tax.rank) %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>%
    psmelt() %>%
    filter(Abundance > abundance.filter)
  #return(abundance.df)
}

GenerateReadSummary <- function(physeq) {
  # Creates a data frames that summarize the total number of reads per RSVs and 
  # per sample (readsPerType). Creates a data frame that just holds the reads
  # per sample (readsPerSample). Creates a summary output of the read 
  # distribution of the samples.
  #
  # Args:
  #   phyloseqObject: An object of class phyloseq.
  #
  # Returns:
  #   List holding readsPerType, readsPerSample, and readDistributionSummary
  
  # readsPerType contains RSVs and Sample names for row names, an nreads column
  # for the total number of reads, a sorted column used as an identifier, and
  # a type column designating if the row is an RSV or a Sample.
  # Total reads per RSV:
  readsPerType <- data.frame(nreads = sort(taxa_sums(physeq),
                                           decreasing = TRUE),
                             sorted = 1:ntaxa(physeq),
                             type = "RSVs")
  # Add total reads per sample:
  readsPerType <- rbind(readsPerType,
                        data.frame(nreads = sort(sample_sums(physeq),
                                                 decreasing = TRUE),
                                   sorted = 1:nsamples(physeq),
                                   type = "Samples"))
  # Create a data frame with just the reads per sample:
  readsPerSample <- data.frame(sum = sample_sums(physeq))
  # Create read distribution summary:
  readDistributionSummary <- summary(readsPerSample)
  
  readSummery <- list("readsPerType" = readsPerType,
                      "readsPerSample" = readsPerSample,
                      "readDistributionSummary" = readDistributionSummary)
}