library(shiny)
dat = readRDS("/home/ajli/shinyApps/carcinogenome/data/MCF10A/data.RDS")
runApp("/home/ajli/shinyApps/carcinogenome/app.R", port = 4848, host = "127.0.0.1")