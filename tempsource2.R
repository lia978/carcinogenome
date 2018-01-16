source("temp.R")


dat = readRDS("./temp/binom.RDS")
runApp(app, port = 3838)