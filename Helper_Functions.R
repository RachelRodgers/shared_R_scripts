# Helper_Functions.Rmd

#----- Colorblind-Friendly Pallets -----#

# Colorblind-friendly pallets for plotting.
cbPaletteGrey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPaletteBlack <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# For fills use
#   scale_fill_manual(values = cbPalletteGrey)

# For line and point colors use
#   scale_colour_manual(values = cbPalleteGrey)

#---------------------------------#

MakeBoxPlot <- function(df, 
                        xVar, 
                        yVar, 
                        label_y = NULL, 
                        label_x = NULL, 
                        statMethod = NULL) {
  # Generate the base plot w/o stat_compare_means
  basePlot <- ggplot(data = df,
                     aes_string(x = xVar,
                                y = yVar)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    ylab(label_y) +
    xlab(label_x) +
    theme_pubr() +
    theme(axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
          axis.text.y = element_text(size = 18))
  # Do we need to add stat_compare_means?
  if (is.null(statMethod)) {
    return(basePlot)
  } else {
    basePlot +
      stat_compare_means(method = statMethod, label.x.npc = 0.5)
  }
}

#---------------------------------#

MakeOrdinationPlot <- function(physeqObj, 
                               ordObj, 
                               colorValues,
                               pointSize = 3.5,
                               shapeValues = NULL,
                               labelColumn = NA,
                               labelSize = 2.5,
                               labelColor = "gray30",
                               shapeVar = NULL, 
                               colorVar = NULL,
                               axesVec = c(1,2),
                               myTitle = NULL,
                               mySubtitle = NULL) {
  
  plot_ordination(physeq = physeqObj,
                  ordination = ordObj,
                  shape = shapeVar,
                  color = colorVar,
                  axes = axesVec) +
    theme_bw() +
    geom_point(size = pointSize) +
    scale_color_manual(values = colorValues) +
    scale_shape_manual(values = shapeValues) +
    ggtitle(myTitle,
            subtitle = mySubtitle) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_text_repel(aes_string(label = labelColumn),
                    size = labelSize,
                    color = labelColor)
}

#---------------------------------#

PlotPhylaPrevalence <- function(prevDF, physeqObj, myTitle = NULL, 
                                mySubtitle = NULL,intercept_y = 0.05, 
                                legendPos = "none") {
  # Calculate total number of samples in the physeqObj
  totalSamples <- nsamples(physeqObj)
  # Add RelativePrev column in prevDF that holds the relative prevalence value
  prevDFMut <- mutate(prevDF, 
                      RelativePrev = prevDF$Prevalence/totalSamples)
  # Generate plot
  ggplot(prevDFMut,
         aes_string(x = "TotalAbundance",
                    y = "RelativePrev",
                    color = "Family")) +
    geom_hline(yintercept = intercept_y, alpha = 0.5,linetype = 2) +
    geom_point(size = 3, alpha = 0.7) +
    scale_x_log10() +
    facet_wrap(~Phylum) +
    xlab("Total Abundance") +
    ylab("Prevalence [Frac. Samples]") +
    ggtitle(myTitle,
            subtitle = mySubtitle) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = legendPos)
}
#---------------------------------#

PlotCommunityComposition <- function(abdDF, taxRank = "Phylum",
                                     facetFormula = NULL,
                                     facetCol = NULL, facetRow = NULL) {
  basePlot <- ggplot(abdDF,
                     aes_string(x = "Sample", y = "Abundance", 
                                fill = taxRank)) +
    geom_bar(stat = "identity", width = 1, color = "grey14") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  # are there facets??
  if (!(is.null(facetFormula))) {
    formula <- as.formula(facetFormula)
    facetPlot <- basePlot +
      facet_wrap(formula, scales = "free", nrow = facetRow, ncol = facetCol)
    return(facetPlot)
  } else {
    return(basePlot)
  }
  
}
#---------------------------------#

PlotAlphaDiversity <- function(df, xVar, yVar, yLabel,
                               statMethod = NULL,
                               alphaPlotTitle = NULL,
                               alphaPlotSubtitle = NULL,
                               facetFormula = NULL,
                               facetCol = NULL,
                               facetRow = NULL) {
  basePlot <- ggplot(df, aes_string(x = xVar, yVar)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    ylab(yLabel) +
    ggtitle(alphaPlotTitle,
            subtitle = alphaPlotSubtitle) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank()) +
    stat_compare_means(method = statMethod, label.x.npc = 0.5, size = 4)
  # are there facets?
  if (!(is.null(facetFormula))) {
    formula <- as.formula(facetFormula)
    facetPlot <- basePlot +
      facet_wrap(formula, ncol = facetCol, nrow = facetRow)
    return(facetPlot)
  } else {
    return(basePlot)
  }
}

#---------------------------------#

# Function to run adonis test on a physeq object and a variable from metadata 
RunAdonis <- function(physeqObj, category, distance) {
  bdist <- phyloseq::distance(physeqObj, distance)
  col <- as(sample_data(physeqObj), "data.frame")[, category]
  # Adonis test
  adonis.bdist <- adonis(bdist ~ col)
  return(adonis.bdist)
}

#---------------------------------#

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

#---------------------------------#

RemoveMissingTaxa <- function(physeq) {
  # Removes any taxa that are missing (have a taxa sum of 0) from a physeq object
  # Error handling
  # Check that a phyloseq object has been passed in
  if ((class(physeq)) != "phyloseq") {
    stop("Please check that you called this function on an object of type 
         \"phyloseq.\"")
  }
  
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
}

#---------------------------------#

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
  
  # readsPerType contains ASVs and Sample names for row names, an nreads column
  # for the total number of reads, a sorted column used as an identifier, and
  # a type column designating if the row is an RSV or a Sample.
  # Total reads per RSV:
  readsPerType <- data.frame(nreads = sort(taxa_sums(physeq),
                                           decreasing = TRUE),
                             sorted = 1:ntaxa(physeq),
                             type = "ASVs")
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
  
  readSummary <- list("readsPerType" = readsPerType,
                      "readsPerSample" = readsPerSample,
                      "readDistributionSummary" = readDistributionSummary)
}

#---------------------------------#

# For community composition plotting.
MakeAbundanceDF <- function(physeq, 
                            taxRank, 
                            abundanceFilter = 0.01,
                            pruneMissing = FALSE) {
  # Creates an abundance data frame at a given taxonomic rank for a phyloseq object
  # for abundance bar plots
  abundance.df <- physeq %>%
    tax_glom(taxrank = taxRank, NArm = pruneMissing) %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>%
    psmelt() %>%
    filter(Abundance > abundanceFilter)
}

#---------------------------------#

TaxRankPrevalence <- function(physeq, taxRank = "Phylum") {
  # Calculate prevalence of taxa at a given taxonomic rank for low prevalence taxon filtering.
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
                                      taxRank,
                                      #"Phylum",
                                      function(df1) {
                                        cbind("Avg_Prevalence" = mean(df1$Prevalence), 
                                              "Total_Prevalence" = sum(df1$Prevalence))
                                      })
  taxRankPrevalence <- list("prevalence_df" = prevalence_df,
                            "prevalence_table" = taxaPrevalence_table)
}

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
  lefseCountsTaxDF <- data.frame(lefseCountsTax)
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
  
  # Merge sample data and maaslin counts by row names
  lefseInput <- merge(lefseSampData, lefseCountsTaxDF, 
                      by = "row.names", all = TRUE)
  lefseInput <- select(lefseInput, -Row.names)
  lefseInput <- t(lefseInput)
  return(lefseInput)
}