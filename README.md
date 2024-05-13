# Hyperion- Imaging Mass Cytometry data analysis

# Installation

To install the stable version from CRAN and bioconductor:

```
install.packages("shiny")
install.packages("tidyverse")
install.packages("data.table")
install.packages("stringr")
install.packages("mclust")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EBImage")
```


Once installed, load the library and run the app:

```
library(shiny)
library(tidyverse)
library(EBImage)
library(stringr)
library(data.table)
library(mclust)

runApp('/imcviewer')
```
