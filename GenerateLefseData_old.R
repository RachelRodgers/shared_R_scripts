
library("phyloseq")
library("tidyverse")

# Function code:
GenerateLefseData <- function(physeqObject, transformCountsToRA = TRUE,
                              categoryColumnName = NULL, sampleColumnName = NULL) {
  # Generates fefse-style input for one category (no sub-categories) from
  # a phyloseq object.
  #
  # Args:
  #   phyloseqObject: An object of class phyloseq.
  #   transformCountsToRA: Specifies whether to transform the counts to relative
  #     abundance (by sample). Default TRUE.
  #   categoryColumnName: The column name of the phyloseq sample data corresponding
  #     to the category that you want to compare counts between in LEfSe.  Function
  #     will stop with an error if this isn't specified.
  #   sampleColumnName: The column name of the phyloseq sample data corresponding 
  #     to the sample names.  If left NULL, will be generated automatically
  #
  # Returns:
  #   A LEfSe-input style data frame that can be written out to file.
  
  # Error handling 
  # Check that the phyloseq library is loaded
  if (!require("phyloseq", character.only = TRUE)) {
    stop("Please make sure the \"phyloseq\" library is loaded.")
  }
  
  # Transform sample counts as relative abundance by sample
  if (transformCountsToRA == TRUE) {
    physeqObject <- transform_sample_counts(physeqObject,
                                            function(OTU) OTU/sum(OTU))
  }
  
  lefseCounts <- data.frame(otu_table(physeqObject))
  lefseCounts <- t(lefseCounts) # samples become rownames
  
  # Add taxonomy to lefseCounts
  lefseCountsTax <- cbind(as.matrix(tax_table(physeqObject)[rownames(lefseCounts), ]), 
                          lefseCounts)
  lefseCountsTaxDF <- data.frame(lefseCountsTax, check.names = FALSE)
  lefseCountsTaxDF <- mutate(lefseCountsTaxDF, 
                             "FullTax" = paste(Kingdom, Phylum, Class,
                                               Order, Family, Genus, Species, 
                                               sep = "|"))
  # Rearrange
  lefseCountsTaxDF <- lefseCountsTaxDF %>%
    dplyr::select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species)) %>%
    dplyr::select(FullTax, everything())
  # Transform back around so we can merge later
  lefseCountsTaxDF <- t(lefseCountsTaxDF)
  colnames(lefseCountsTaxDF) <- lefseCountsTaxDF[1, ]
  lefseCountsTaxDF <- lefseCountsTaxDF[-1, ]
  
  # Get the sample data with the metadata variables of choice
  lefseSampData <- data.frame(sample_data(physeqObject))
  # Modify the sample data to contain sample names in a column and category of interest
  # First check if we have a category column name, stop otherwise 
  if (is.null(categoryColumnName)) {
    stop("Need to specify a category")
  }
  
  # Next check for a sample column name, if not, create one
  if (is.null(sampleColumnName)) {
    lefseSampData <- lefseSampData %>%
      rownames_to_column(var = "Sample") %>%
      mutate("SampleID" = Sample) %>%
      column_to_rownames(var = "Sample")
    
    lefseSampData <- lefseSampData[ , c("SampleID", categoryColumnName)]
  } else {
    lefseSampData <- lefseSampData[, c(sampleColumnName, categoryColumnName)]
  }
  
  # Merge sample data and lefse counts by row names
  lefseInput <- merge(lefseSampData, lefseCountsTaxDF, 
                      by = "row.names", all = TRUE)
  lefseInput <- select(lefseInput, -Row.names)
  lefseInput <- t(lefseInput)
  return(lefseInput)
}

# Test data
#myPhyseq <- readRDS("./physeqEpg5KOvWT.RDS")

# Check which columns we need for the lefse data
#sample_variables(myPhyseq)

# We will need the SampleName column (contains the sample identifiers) and the
#   Genotype column which contains our variables of interest (tells us whether
#   a sample is of genotype KO or WT).

# Note: lefse can normalize the counts on its own, so you could set the transformCountsToRA parameter to FALSE

#lefseData <- GenerateLefseData(physeqObject = myPhyseq,
                               #transformCountsToRA = TRUE,
                               #categoryColumnName = "Genotype", 
                               #sampleColumnName = "SampleName")

# Save the lefse data
#write.table(x = lefseData, file = "lefseData_RelAbd.txt", 
            #quote = FALSE, sep = "\t", col.names = FALSE)
