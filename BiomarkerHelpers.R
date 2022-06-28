# BiomarkerHelpers.R

#---------------------------------#
# DESeq2 - specific functions (old):
GetBiomarkers <- function(physeq, 
                          groupVar,
                          numerator, 
                          denominator,
                          alpha = 0.05) { 
  # Create the DESeq significance table used for plotting.
  # Returns a list of two objects - the DESeq2 results and the Significance Table (sigTable)
  # Each component can be accessed with dollar sign ($) notation
  
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
  biomarkerResults <- list("results" = ddsResults, "biomarkerTable" = ddsSignificantResults)
  #return(ddsSignificantResults)
  #biomarkerResults
}

# Horizontal Biomarker Plot

CreateBiomarkerPlot <- function(biomarkerTable) {
  biomarkerPlot <- ggplot(data.frame(biomarkerTable),
                          aes(x = Family,
                              y = log2FoldChange,
                              color = Genus)) +
    geom_point(aes(size = log10(biomarkerTable$baseMean)), 
               alpha = 0.7) +
    geom_hline(yintercept = 0, lwd = 1.5) +
    ggtitle("default plot title") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 315, hjust = 0)) +
    labs(size = "log10(base mean)")
  return(biomarkerPlot)
}

# Volcano Biomarker Plot

CreateBiomarkerVolcano <- function(physeq, biomarkerResults, alpha = 0.05) {
  # Set up
  taxonomyTable <- data.table(data.frame(as(tax_table(physeq), "matrix")), 
                              keep.rownames = TRUE)
  setnames(taxonomyTable, "rn", "OTU")
  setkeyv(taxonomyTable, "OTU")
  
  resultsDataTable <- data.table(as(biomarkerResults$results, "data.frame"),
                                 keep.rownames = TRUE)
  setnames(resultsDataTable, "rn", "OTU")
  setkeyv(resultsDataTable, "OTU")
  
  resultsDataTable <- taxonomyTable[resultsDataTable]
  resultsDataTable <- resultsDataTable %>%
    filter(., padj != "NA") %>%
    mutate(., Significant = padj <alpha)
  
  # Create volcano plot object
  volcano <- ggplot(resultsDataTable,
                    aes(x = log2FoldChange,
                        y = -log10(padj))) +
    geom_point(data = subset(resultsDataTable, resultsDataTable$Significant == FALSE), color = "grey") +
    geom_point(data = subset(resultsDataTable, resultsDataTable$Significant == TRUE), 
               aes(color = Phylum, size = baseMean)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_hline(yintercept = -log10(alpha)) +
    ggtitle("default plot title") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12)) 
  return(volcano)
}

## TEST ##
# how to use the DESeq2 specific functions
#
# Physeq object: physeqWTDay21
#
# Get biomarker results and significance table from DESeq2 comparing the
#   Treatment groups with No_Antibiotics as the numerator (top_), Antibiotics 
#   as the denominator (bottom) of the Horizontal Biomarkers plot:
# 
# biomarkers_WTDay21_NoAbx_Abx <- GetBiomarkers(physeqWTDay21, 
#                                              "Treatment", 
#                                              "No_Antibiotics", "Antibiotics")
# Create Horizontal plot
# biomarkerPlot_WTDay21_NoAbx_Abx <- CreateBiomarkerPlot(biomarkers_WTDay21_NoAbx_Abx$biomarkerTable)
# biomarkerPlot_WTDay21_NoAbx_Abx
# 
# Create volcano plot
# biomarkerVolcano_WTDay21_NoAbx_Abx <- CreateBiomarkerVolcano(physeqWTDay21,
#                                                             biomarkers_WTDay21_NoAbx_Abx)
# biomarkerVolcano_WTDay21_NoAbx_Abx

#---------------------------------#

# Functions to quickly make interactive volcano plots:

GenerateDESeqResults <- function(physeq, variable, numerator, denominator) {
  # Returns DESeq Results as Formal Class "DESeqResults"
  # Create formula from string variable
  formula <- as.formula(paste("~", variable, sep = " "))
  # Convert to deseq data set object
  dds <- phyloseq_to_deseq2(physeq, design = formula)
  # Run analysis
  ddsAnalysis <- DESeq(dds, test = "Wald", fitType = "local", betaPrior = FALSE)
  # Extract and format results
  ddsResults <- results(ddsAnalysis,
                        contrast = c(variable, numerator, denominator)) 
}

#mcols(ddsResults)

GenerateDESeqResultsTable <- function(physeq, ddsResults, sigThreshold = 0.05) {
  # Returns data frame
  # From the DESeq results generated by GenerateDESeqResults, create a
  #   results data table that includes the taxonomy information and a column
  #   indicating whether results for each taxon are significant.
  
  # Extract taxonomy table:
  taxTable <- data.table(data.frame(as(tax_table(physeq), "matrix")),
                         keep.rownames = TRUE)
  setnames(taxTable, "rn", "OTU")
  setkeyv(taxTable, "OTU")
  
  # Extract DESeq results as a data frame:
  resDT <- data.table(as(ddsResults, "data.frame"),
                      keep.rownames = TRUE)
  
  setnames(resDT, "rn", "OTU")
  setkeyv(resDT, "OTU")
  
  # Combine taxonomy information with the results table:
  resDT <- taxTable[resDT]
  resDT <- resDT %>%
    filter(padj != "NA") %>%
    mutate(Significant = padj < sigThreshold)
}

PlotStaticVolcano <- function(physeq,
                              resultsDataTable,
                              sigThreshold, # usually 0.05 (match sigThreshold in GenerateDESeqResultsTable())
                              plotTitle = NULL) {
  # Returns ggplot object from results data frame generated from
  #   GenerateDESeqResultsTable()
  
  # Create volcano plot object
  volcano <- ggplot(resultsDataTable,
                    aes(x = log2FoldChange,
                        y = -log10(padj),
                        label1 = Family,
                        label2 = Genus,
                        label3 = Species)) +
    geom_point(data = subset(resultsDataTable,
                             resultsDataTable$Significant == FALSE),
               color = "grey") +
    geom_point(data = subset(resultsDataTable,
                             resultsDataTable$Significant == TRUE),
               aes(color = Phylum, size = baseMean)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_hline(yintercept = -log10(sigThreshold)) +
    ggtitle(plotTitle) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12))
}

## Example Usage ##

#resultsDESeq <- GenerateDESeqResults(physeq = physeqBacteria,
#                                     variable = "Condition",
#                                     numerator = "HHC",
#                                     denominator = "CRPS")


#resTable <- GenerateDESeqResultsTable(physeq = physeqBacteria,
#                                      ddsResults = resultsDESeq)


#volcano <- PlotStaticVolcano(physeq = physeqBacteria,
#                             resultsDataTable = resTable,
#                             plotTitle = "Differentially Abundant Taxa \nHHC Relative to CRPS")

#ggplotly(volcano, tooltip = c("Phylum", "Genus", "Species",
#                              "log2FoldChange", "baseMean"))

#---------------------------------#
