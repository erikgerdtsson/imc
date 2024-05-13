# imc
Hyperion- Imaging Mass Cytometry data analysis
Installation
To install the stable version from CRAN:

install.packages("shiny")
Getting Started
Once installed, load the library and run an example:

library(shiny)
# Launches an app, with the app's source code included
runExample("06_tabsets")
# Lists more prepackaged examples
runExample()
For more examples and inspiration, check out the Shiny User Gallery.


library(shiny)
library(tidyverse)
library(EBImage)
library(stringr)
library(data.table)
library(mclust)
