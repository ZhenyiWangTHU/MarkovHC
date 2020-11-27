# Install requirements: parallel, doParallel, dbscan, igraph, Matrix, Rcpp, plyr, dplyr, doBy, ggraph
# Require R version >= 3.5

if (as.numeric(R.version$major) <= 3 && as.numeric(R.version$minor) < 5) {
  # R < 3.5
  message("Require R version >= 3.5")
} else {
  # R >= 3.5
  message("Installing required packages from cran.")
  install.packages(c("parallel", "doParallel", "dbscan", "igraph", "Matrix", "Rcpp", "plyr", "dplyr", "doBy", "ggraph"))
}

# Check that packages installation went smoothly.
if (!requireNamespace("parallel", quietly = TRUE)) {stop("Failed to install required package 'parallel' from cran.")}
if (!requireNamespace("doParallel", quietly = TRUE)) {stop("Failed to install required package 'doParallel' from cran.")}
if (!requireNamespace("dbscan", quietly = TRUE)) {stop("Failed to install required package 'dbscan' from cran.")}
if (!requireNamespace("igraph", quietly = TRUE)) {stop("Failed to install required package 'igraph' from cran.")}
if (!requireNamespace("Matrix", quietly = TRUE)) {stop("Failed to install required package 'Matrix' from cran.")}
if (!requireNamespace("Rcpp", quietly = TRUE)) {stop("Failed to install required package 'Rcpp' from cran.")}
if (!requireNamespace("plyr", quietly = TRUE)) {stop("Failed to install required package 'plyr' from cran.")}
if (!requireNamespace("dplyr", quietly = TRUE)) {stop("Failed to install required package 'dplyr' from cran.")}
if (!requireNamespace("doBy", quietly = TRUE)) {stop("Failed to install required package 'doBy' from cran.")}
if (!requireNamespace("ggraph", quietly = TRUE)) {stop("Failed to install required package 'ggraph' from cran.")}

# Install MarkovHC (and devtools, if necessary)
if (requireNamespace("devtools", quietly = TRUE)) {
  message("Installing MarkovHC")
  devtools::install_github(repo="ZhenyiWangTHU/MarkovHC")
} else {
  message("Installing devtools")
  install.packages("devtools")
  message("Installing MarkovHC")
  devtools::install_github(repo="ZhenyiWangTHU/MarkovHC")
}

# Check that MarkovHC installed.
if (requireNamespace("MarkovHC", quietly = TRUE)) {
  message("MarkovHC installed successfully!")
  message('You can load it by typing: library("MarkovHC")')
  message('Try "?MarkovHC" for starting tips.')
} else {
  message("Something went wrong. It doesn't seem that MarkovHC installed correctly.")
}
