# R_helper_functions.R
# From Eli Roberson

##########################
# Package load functions #
##########################

# Note: require() returns logical - 
# FALSE if the requested package is not found and TRUE if the package is loaded

load_bioc_libraries = function( libName )
{
  if ( !require( libName, character.only = TRUE ) )
  {
    source( "http://www.bioconductor.org/biocLite.R" )
    biocLite( libName, suppressUpdates=TRUE )
    
    if ( !require( libName, character.only = TRUE ) )
    {
      stop( "Couldn't install ", libName, " from Bioconductor. Please install manually." )
    }
  }
}


# This function will try CRAN first.
load_R_libraries = function( libName )
{
  if ( !require( libName, character.only=TRUE ) )
  {
    install.packages( libName, repos="http://cran.wustl.edu" )
    
    if ( !require( libName, character.only=TRUE ) )
    {
      stop( "Couldn't install ", libName, " from R packages. Please install manually." )
    }
    
  }
}


