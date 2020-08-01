# Installing MarkovHC:

## Quick Installation

### 1. Install R and RStudio

**MarkovHC** is a package written in R and designed to be used in the RStudio interactive environment.

R can be obtained and installed from [https://cran.rstudio.com](https://cran.rstudio.com). 

MarkovHC has only been tested under R version 3.5.1 and 3.6.1, and may fail to install under R 3.4.x until some of its dependencies are updated.

Following installation of R, Rstudio can be obtained and installed from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/). The free version of RStudio Desktop is sufficient.

### 2. Install required R packages and MarkovHC

We wrote a script that will attempt to install all requirements for MarkovHC and then MarkovHC itself. Start RStudio and then run the installation script from the console:

```source("https://raw.githubusercontent.com/farrellja/URD/master/URD-Install.R")```

## Troubleshooting

### udunits

URD's dependencies require the installation of the **udunits2** package, which depends on the udunits2 library. This is part of the Windows and Mac OS X operating systems, but often requires manual installation on Linux distributions:

```
sudo apt-get install libudunits2-dev
```

### DLL error

R has limit on the number of DLLs that can be loaded by linked packages. If you receive a **maximal number of DLLs reached** error during installation, then you can increase the number of simultaneously allowed DLLs from 100 to 200 by modifying the *.Renviron* file.

On Mac, from the terminal, run:
```echo "R_MAX_NUM_DLLS=200" >> ~/.Renviron```

### Manual installation

If something went wrong previously, you may have try installing some of URD's dependencies manually:

###### A. Install and attach the *devtools* package

*devtools* is required to compile and install URD, since it's distributed through GitHub.

```install.packages("devtools")```
     
###### B. Install required Bioconductor packages

Because these packages are deposited in Bioconductor and not CRAN, they must be installed manually before installing URD (otherwise its installation will fail.)

```source("https://bioconductor.org/biocLite.R")```
```biocLite(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny'))```

Optionally, additional packages that are required for specific functions can be installed.

```biocLite(c('sva', 'rhdf5', 'scran'))```

If installing on a cluster, where you do not have write permissions to the main R libraries (for instance, if you receive errors like "ERROR: no permission to install to directory", you may need to specify a library location, such as:

```biocLite(c('sva', 'rhdf5', 'scran'), lib=[PATH TO DIRECTORY TO STORE R LIBRARIES WHERE YOU HAVE PERMISSION TO WRITE])```
     
###### C. Install URD

*URD* can be installed directly from the GitHub repository

```library(devtools)```
```install_github(repo="farrellja/URD")```
