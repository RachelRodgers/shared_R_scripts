# DESeq2_Helper_Functions.R


GetBiomarkersTable <- function(physeq, 
                               groupVar,
                               numerator, 
                               denominator,
                               alpha = 0.05) { 
  # Create formula from string variable
  formula <- as.formula(paste("~", groupVar, sep = " "))
  # Convert physeq object to DESeq Data Set object
  dds <- phyloseq_to_deseq2(physeq = physeq, design = formula)
  
  # Run DESeq analysis
  ddsAnalysis <- DESeq(dds, test = "Wald", fitType = "local", betaPrior = FALSE)
  # Extract Results
  ddsResults <- results(ddsAnalysis,
                        contrast = c(groupVar, numerator, denominator),
                        cooksCutoff = FALSE)
  # Create table of significant results
  ddsSignificantResults <- ddsResults[which(ddsResults$padj < alpha), ]
  # Add taxonomy to the table
  ddsSignificantResults <- cbind(as(ddsSignificantResults, "data.frame"),
                                 as(tax_table(physeq)[rownames(ddsSignificantResults), ], "matrix"))
  # Get the maximum log2FC for each Family, then sort in decreasing order
  maxFC <- tapply(ddsSignificantResults$log2FoldChange,
                  ddsSignificantResults$Family,
                  function(x) max(x))
  maxFC <- sort(maxFC, decreasing = TRUE)
  # Change significance table Species column to Families, factorize the column and assign levels by decreasing max log2FC
  ddsSignificantResults$Species <- factor(as.character(ddsSignificantResults$Family),
                                          levels = names(maxFC))
  return(ddsSignificantResults)
}


