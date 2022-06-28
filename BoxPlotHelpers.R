# BoxPlotHelpers.R

#----- Box Plots -----##

MakeBoxPlot <- function(df, xVar, yVar, label_y = NULL, label_x = NULL, 
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


PrettyBox <- function(df, xVar, yVar, statMethod, plotTitle,
                      colorVals = NULL,
                      label_x = NULL, label_y = NULL, legendPos = "none",
                      facet_formula = NULL, facet_rows = NULL, facet_cols = NULL) {
  
  basePlot <- ggplot(data = df,
                     aes_string(x = xVar, y = yVar, fill = xVar)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5) +
    scale_fill_manual(values = colorVals) +
    labs(x = label_x, y = label_y, title = plotTitle) +
    theme_pubr() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = 0.5),
          strip.text = element_text(size = 12),
          legend.title = element_blank(), 
          legend.position = legendPos,
          legend.text = element_text(size = 12)) +
    stat_compare_means(method = statMethod, label.x.npc = 0)
  
  if (!is.null(facet_formula)) {
    
    facetPlot <- basePlot +
      facet_wrap(as.formula(facet_formula), nrow = facet_rows, ncol = facet_cols,
                 scales = "free")
    
    return(facetPlot)
    
  } else {
    return(basePlot)
  }
  
}