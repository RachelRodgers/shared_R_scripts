# Helper_Functions.Rmd

#----- Colorblind-Friendly Pallets -----#

# Colorblind-friendly pallets for plotting.
cbPaletteGrey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPaletteBlack <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# For fills use
#   scale_fill_manual(values = cbPalletteGrey)

# For line and point colors use
#   scale_colour_manual(values = cbPalleteGrey)

#----- Statistical Functions -----#

# Function to run adonis test on a physeq object and a variable from metadata 
RunAdonis <- function(physeqObj, category, distance) {
  bdist <- phyloseq::distance(physeqObj, distance)
  col <- as(sample_data(physeqObj), "data.frame")[, category]
  # Adonis test
  adonis.bdist <- adonis(bdist ~ col)
  return(adonis.bdist)
}

#----- Common Tasks Code -----#

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

#----- Common QC Tasks -----#

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

#----- Common 16S Analysis Tasks -----#

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


#----- Ordination Tasks -----#

EasyOrd <- function(physeqObj, ordObj, colorVals, pointSize = 1.5,
                    shapeVar = NULL, colorVar = NULL, myTitle = NULL,
                    mySubtitle = NULL, myLegendTitle = NULL) {
  plot_ordination(physeq = physeqObj,
                  ordination = ordObj,
                  shape = shapeVar,
                  color = colorVar) +
    scale_color_manual(values = colorVals) +
    labs(color = myLegendTitle) +
    theme_bw() +
    geom_point(size = pointSize) +
    ggtitle(myTitle, subtitle = mySubtitle) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12)) +
    stat_ellipse(type = "norm")
}


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