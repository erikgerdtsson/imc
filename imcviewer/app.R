#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# update the project directory before running.
data.directory <- "/my_path_here/"




library(shiny)
library(tidyverse)
library(EBImage)
library(stringr)
library(data.table)
library(mclust)

ggname <- function(x) {
  if (class(x) != "character") {
    return(x)
  }
  y <- sapply(x, function(s) {
    if (!grepl("^`", s)) {
      s <- paste("`", s, sep = "", collapse = "")
    }
    if (!grepl("`$", s)) {
      s <- paste(s, "`", sep = "", collapse = "")
    }
  })
  y
}


cleanChannelName <- function(x) {
  string.to.remove <- str_split(x, "_", simplify = T)[, 1]
  channelName <- gsub(paste0(string.to.remove, "_"), "", x)
  return(channelName)
}



intensity.files <- list.files(
  path = paste0(data.directory, "intensities/"), pattern = ".csv",
  all.files = TRUE, full.names = F, no.. = TRUE
)

regionprops.files <- list.files(
  path = paste0(data.directory, "regionprops/"), pattern = ".csv",
  all.files = TRUE, full.names = F, no.. = TRUE
)



data_list.intensities <- lapply(intensity.files, function(x) data.table(fread(paste0(data.directory, "intensities/", x))) %>% mutate(name = gsub(".csv", "", x))) %>% rbindlist(., fill = TRUE)
data_list.regionprops <- lapply(regionprops.files, function(x) data.table(fread(paste0(data.directory, "regionprops/", x))) %>% mutate(name = gsub(".csv", "", x))) %>% rbindlist(., fill = TRUE)

data <- left_join(data_list.intensities, data_list.regionprops)


rois <- list.dirs(
  path = paste0(data.directory, "histocat/"),
  full.names = F, recursive = F
)

channels <- list.files(
  path = paste0(data.directory, "histocat/", rois[1]), pattern = ".tiff",
  all.files = TRUE, full.names = F, no.. = TRUE
) %>% gsub(".tiff", "", .)

if (grep("mask", channels) > 0) {
  channels <- channels[-grep("mask", channels)]
}


names(channels) <- lapply(channels, cleanChannelName)
channels <- channels[order(names(channels))]


ui <- fluidPage(

  # Application title
  titlePanel("IMC display"),

  # Sidebar with a select input for the image
  sidebarLayout(
    sidebarPanel(
      selectInput("roi", "ROI:", rois),
      selectInput("green", "Green channel:", channels),
      numericInput("gain_green", "max", 10, min = 0, max = 100),
      checkboxInput("mask", "show mask", FALSE),
      numericInput("threshold", "threshold", 10, min = 0, max = 100),
      actionButton("go", "Auto Threshold"),
      textOutput("text"),
      selectInput("red", "Red channel:", channels),
      numericInput("gain_red", "max", 10, min = 0, max = 100),
      selectInput("blue", "Blue channel:", channels),
      numericInput("gain_blue", "max", 10, min = 0, max = 100),
      plotOutput("hist")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Green with mask", displayOutput("widget2", height = 800)),
        tabPanel("Composite", displayOutput("widget", height = 800)),
        tabPanel("plot", checkboxInput("arcsin", "archsin transform", FALSE), plotOutput("plot", height = 1000), )
      )
    )
  )
)

server <- function(input, output, session) {
  img <- reactive({
    path_roi <- paste0(data.directory, "histocat/", input$roi, "/")
    green <- readImage(paste0(path_roi, input$green, ".tiff")) / (input$gain_green)
    red <- readImage(paste0(path_roi, input$red, ".tiff")) / (input$gain_red)
    blue <- readImage(paste0(path_roi, input$blue, ".tiff")) / (input$gain_blue)
    comp <- rgbImage(green = green, red = red, blue = blue)

    if (input$mask) {
      mask <- readImage(paste0(path_roi, input$roi, "_mask.tiff"))
      new <- paintObjects(mask, comp, col = "purple")
    } else {
      comp
    }
  })

  estimatedTreshold <- eventReactive(input$go, {
    green <- cleanChannelName(input$green)
    df <- data %>%
      subset(name == input$roi) %>%
      select(., green) %>%
      as.matrix()
    cur_model <- Mclust(asinh(df), G = 2)
    max1 <- max(df[cur_model$classification == 1])
    max2 <- max(df[cur_model$classification == 2])
    auto_thresh <- ifelse(max1 < max2, max1, max2) * 2
    updateNumericInput(session, "threshold", value = auto_thresh)
    return(auto_thresh)
  })

  img2 <- reactive({
    path_roi <- paste0(data.directory, "histocat/", input$roi, "/")
    green <- readImage(paste0(path_roi, input$green, ".tiff")) / (input$gain_green)
    comp <- rgbImage(green = green)


    if (input$mask) {
      mask <- readImage(paste0(path_roi, input$roi, "_mask.tiff"))
      mask <- round(mask * 2^16)
      green.channel <- cleanChannelName(input$green)
      t <- data %>%
        filter(name == input$roi) %>%
        filter(!!sym(green.channel) > input$threshold)
      mask[!mask %in% t$Object] <- 0
      paintObjects(mask, comp, col = "white")
    } else {
      comp
    }
  })

  max_x <- reactive({
    red <- cleanChannelName(input$red)
    df <- data %>%
      select(., red) %>%
      as.matrix()
    quantile(df, 0.99)
  })

  max_y <- reactive({
    green <- cleanChannelName(input$green)
    df <- data %>%
      select(., green) %>%
      as.matrix()
    quantile(df, 0.99)
  })

  output$widget <- renderDisplay({
    display(img())
  })

  output$text <- renderText({
    green.channel <- cleanChannelName(input$green)
    t <- data %>%
      filter(name == input$roi) %>%
      filter(!!sym(green.channel) > input$threshold)
    t2 <- data %>%
      filter(name == input$roi) %>%
      filter(!!sym(green.channel) <= input$threshold)
    n <- nrow(t)
    snr <- mean(t[[green.channel]]) / mean(t2[[green.channel]])
    paste("n=", n, " ,SNR=", round(snr))
  })

  output$estimatedTreshold <- renderText({
    round(estimatedTreshold(), 1)
  })

  output$widget2 <- renderDisplay({
    display(img2())
  }, )


  output$hist <- renderPlot({
    green <- cleanChannelName(input$green)
    df <- data %>%
      subset(name == input$roi) %>%
      select(., green) %>%
      as.matrix()

    d <- density(asinh(df))
    plot(d, main = "signal int.")

    abline(v = asinh(input$threshold), col = "red", lwd = 3)
    abline(v = asinh(round(estimatedTreshold(), 1)), col = "green", lwd = 2)
  })





  output$plot <- renderPlot({
    if (!input$arcsin) {
      red <- ggname(cleanChannelName(input$red))
      green <- ggname(cleanChannelName(input$green))
      str(red)
      maxx <- max_x()
      maxy <- max_y()


      ggplot(data, aes_string(x = red, y = green)) +
        geom_point(alpha = 0.1, size = 0.1, color = "red", na.rm = T) +
        facet_wrap(~name) +
        xlab(red) +
        ylab(green) +
        geom_hline(aes(yintercept = input$threshold), color = "blue", linetype = "dotted", size = 0.5) +
        xlim(0, maxx) +
        ylim(0, maxy) +
        theme(aspect.ratio = 1)
    } else {
      xcolumn <- cleanChannelName(input$red)
      ycolumn <- cleanChannelName(input$green)

      if (xcolumn != ycolumn) {
        ggplot(data %>% select(c(xcolumn, ycolumn, "name")) %>% set_names(c("x", "y", "name")), aes(asinh(x / 1), asinh(y / 1))) +
          geom_point(alpha = 0.1, size = 0.1, color = "red") +
          facet_wrap(~name) +
          geom_hline(aes(yintercept = asinh(input$threshold / 1)), color = "blue", linetype = "dotted", size = 0.5) +
          xlab(xcolumn) +
          ylab(ycolumn) +
          theme(aspect.ratio = 1)
      }
    }
  })
}


# Run the application
shinyApp(ui = ui, server = server)
