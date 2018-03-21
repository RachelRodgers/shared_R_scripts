# function_playground.R

library(dplyr)
library(tibble)

test_table_taxaAreRows <- read.delim("test_otutable.txt",
                                     header = TRUE,
                                     sep = "\t")
row.names(test_table_taxaAreRows) <- test_table_taxaAreRows$X
test_table_taxaAreRows <- test_table_taxaAreRows[-1]

test_table_taxaAreColumns <- data.frame(t(test_table_taxaAreRows))

buildRareCurves <- function(physeq, taxa_are_rows) {

  otuTable <- data.frame(otu_table(physeq))

  # Transpose otuTable if taxa_are_rows(physeq) evals to FALSE
  if (!taxa_are_rows(physeq)) {
    otuTable <- data.frame(t(otuTable))
  }
  # taxa/OTUs are the row names
  # Change OTUs from row names to its own column called OTU
  otuTable <- rownames_to_column(otuTable, var = "OTU")
  
  # Vectors to hold loop output
  sampleVector <- vector(mode = "numeric")
  nonuniqueTaxa <- vector(mode = "numeric")
  uniqueTaxa <- vector(mode = "numeric")
  
  # Fill the vectors
  for (j in 2:ncol(otuTable)) {
    currentSampleIndex <- j - 1 # since j starts at 2, not 1
    sampleVector[currentSampleIndex] <- currentSampleIndex
    
    filteredTable <- filter(otuTable, otuTable[, j] > 0) # Remove missing OTUs
    nonuniqueTaxa <- c(nonuniqueTaxa, filteredTable$OTU) # Put remaining OTUs in vector
    uniqueTaxa[currentSampleIndex] <- length(unique(nonuniqueTaxa)) # Count the unique taxa
  }
  
  rareCurveDF <- data.frame(cbind("Sample" = sampleVector,
                                  "UniqueTaxa" = uniqueTaxa))
}

test_taxaAreRows <- buildRareCurves(test_table_taxaAreRows, taxa_are_rows = TRUE)
test_taxaAreColumns <- buildRareCurves(test_table_taxaAreColumns, taxa_are_rows = FALSE)
