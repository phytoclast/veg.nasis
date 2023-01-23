library(soilDB)
library(stringr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

veg.plot <- soilDB::get_vegplot_from_NASIS_db()

veg.data <- soilDB::get_veg_data_from_NASIS_db()

veg.loc <- soilDB::get_vegplot_location_from_NASIS_db()

veg.spp <- soilDB::get_vegplot_species_from_NASIS_db()

veg.data.veg <- veg.data$veg

veg.data.trans <- veg.data$vegtransect


ht.metric <- function(ft){
  round(ft*0.3048,1)
}
ht.medieval <- function(m){
  round(m/0.3048,1)
}
diam.metric <- function(inch){
  round(inch*2.54,0)
}
diam.medieval <- function(cm){
  round(cm/2.54,1)
}


syns <- read.csv('data/plants/m.ac.csv')