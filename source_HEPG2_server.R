library(shiny)
dat = readRDS("/home/ajli/shinyApps/carcinogenome/data/HEPG2/data.RDS")
runApp("/home/ajli/shinyApps/carcinogenome/app.R", port = 3838, host = "127.0.0.1")