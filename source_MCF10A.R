dat = readRDS("./data/MCF10A/data.RDS")
source("app.R")
runApp(app, port = 4848)