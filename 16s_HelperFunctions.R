#16S_HelperFunctions.R

################################################################################
# Colorblind-friendly pallets for plotting.
cbPaletteGrey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPaletteBlack <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# For fills use
#   scale_fill_manual(values = cbPalletteGrey)

# For line and point colors use
#   scale_colour_manual(values = cbPalleteGrey)



################################################################################
# Helper functions.

## Functions for basic manipulation of phyloseq object. ##

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


## Functions for generating data frames from phyloseq objects ##

# For pre-processing based on read distributions.
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

# For community composition plotting.
MakeAbundanceDF <- function(physeq, taxRank, abundanceFilter = 0.01) {
  # Creates an abundance data frame at a given taxonomic rank for a phyloseq object
  # for abundance bar plots
  abundance.df <- physeq %>%
    tax_glom(taxrank = taxRank) %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>%
    psmelt() %>%
    filter(Abundance > abundanceFilter)
  #return(abundance.df)
}

TaxRankPrevalence <- function(physeq, taxRank = "Phylum") {
  # Create a named vector where each element name is an OTU sequeunce, and each value
  #   is the number of sequences in which that OTU is present (max value is the total
  #   number of sequences).
  prevalence_vector <- apply(X = otu_table(physeq),
                             MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no =2),
                             FUN = function(x) {sum(x > 0)})
  # Generate a prevalence dataframe that also adds a TotalAbundance column (the total
  #   number of reads for that OTU across all samples) and the taxonomy information 
  #   for each OTU.
  prevalence_df <- data.frame(Prevalence = prevalence_vector,
                              TotalAbundance = taxa_sums(physeq),
                              tax_table(physeq))
  # Create a new data frame that displays, for a given taxonomic rank, the average
  #   number of samples in which that taxon is present, and the total number of samples
  #   in which that taxon is present.
  taxaPrevalence_table <- plyr::ddply(prevalence_df,
                                      "Phylum",
                                      function(df1) {
                                        cbind("Avg_Prevalence" = mean(df1$Prevalence), 
                                              "Total_Prevalence" = sum(df1$Prevalence))
                                      })
  taxRankPrevalence <- list("prevalence_df" = prevalence_df,
                            "prevalence_table" = taxaPrevalence_table)
}

## Functions for creating of plotting data frames. ##

buildRareCurves <- function(physeq) {
  # Make a data frame from the OTU table
  otu_table <- data.frame(otu_table(physeq))
  # Change OTUs from row names to a column
  otu_table <- rownames_to_column(otu_table, var = "OTU")
  
  # Vectors to hold for loop output
  sample_vector <- vector(mode = "numeric")
  nonunique_taxa <- vector(mode = "numeric")
  unique_taxa <- vector(mode = "numeric")
  
  # Fill vectors
  for (i in 2:ncol(otu_table)) {
    current_sample_index <- i - 1
    sample_vector[current_sample_index] <- current_sample_index
    
    filtered_df <- filter(otu_table, otu_table[, i] > 0)
    nonunique_taxa <- c(nonunique_taxa, filtered_df$OTU)
    unique_taxa[current_sample_index] <- length(unique(nonunique_taxa))
  }
  
  # Bind vectors into a data frame (this is the return value)
  rare_curve_df <- data.frame(cbind("Sample" = sample_vector,
                                    "UniqueTaxa" = unique_taxa))
  
}


# Calculate prevalence of taxa at a given taxonomic rank for low prevalence taxon filtering. #
TaxRankPrevalence <- function(physeq, taxRank = "Phylum") {
  # Create a named vector where each element name is an OTU sequeunce, and each value
  #   is the number of sequences in which that OTU is present (max value is the total
  #   number of sequences).
  prevalence_vector <- apply(X = otu_table(physeq),
                             MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no =2),
                             FUN = function(x) {sum(x > 0)})
  # Generate a prevalence dataframe that also adds a TotalAbundance column (the total
  #   number of reads for that OTU across all samples) and the taxonomy information 
  #   for each OTU.
  prevalence_df <- data.frame(Prevalence = prevalence_vector,
                              TotalAbundance = taxa_sums(physeq),
                              tax_table(physeq))
  # Create a new data frame that displays, for a given taxonomic rank, the average
  #   number of samples in which that taxon is present, and the total number of samples
  #   in which that taxon is present.
  taxaPrevalence_table <- plyr::ddply(prevalence_df,
                                      "Phylum",
                                      function(df1) {
                                        cbind("Avg_Prevalence" = mean(df1$Prevalence), 
                                              "Total_Prevalence" = sum(df1$Prevalence))
                                      })
  taxRankPrevalence <- list("prevalence_df" = prevalence_df,
                            "prevalence_table" = taxaPrevalence_table)
}