# Installing MarkovHC:

## Quick Installation

### 1. Install R and RStudio

**MarkovHC** is a package written in R and designed to be used in the RStudio interactive environment.

R can be obtained and installed from [https://cran.rstudio.com](https://cran.rstudio.com). 

MarkovHC has only been tested under R version 3.5.1 and 3.6.1, and may fail to install under R 3.4.x until some of its dependencies are updated.

Following installation of R, Rstudio can be obtained and installed from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/). The free version of RStudio Desktop is sufficient.

### 2. Install required R packages and MarkovHC

We wrote a script that will attempt to install all requirements for MarkovHC and then MarkovHC itself. Start RStudio and then run the installation script from the console:

```source("https://github.com/ZhenyiWangTHU/MarkovHC/blob/master/MarkovHC-install.R")```

### Manual installation

If something went wrong previously, you may have try installing some of MarkovHC's dependencies manually:

###### A. Install and attach the *devtools* package

*devtools* is required to compile and install MarkovHC, since it's distributed through GitHub.

```install.packages("devtools")```
     
###### B. Install required packages

Because these packages must be installed before installing MarkovHC (otherwise its installation will fail.)

``` install.packages(c("parallel", "doParallel", "dbscan", "igraph", "Matrix", "Rcpp", "plyr", "dplyr", "doBy", "ggraph"))```
     
###### C. Install MarkovHC

*MarkovHC* can be installed directly from the GitHub repository

```library(devtools)```  
```install_github(repo="ZhenyiWangTHU/MarkovHC",subdir = "/MarkovHC")```

Or download MarkovHC_1.0.0.tar.gz and install in R

```install.packages('MarkovHC_1.0.0.tar.gz')```
