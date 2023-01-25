library(soilDB)
library(stringr)
library(dplyr)

# remotes::install_github("natearoe/ecositer", dependencies = FALSE)
# library(ecositer)
# test <- veg_summary(veg_df = ecositer::vegetation_dataframe)
# 
# test.R022AB006CA <- test$R022AB006CA$Raw_data



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

veg.plot <- readRDS('data/veg.plot.RDS')
veg.data <-  readRDS('data/veg.data.RDS')
veg.loc <-  readRDS('data/veg.loc.RDS')
veg.spp <-  readRDS('data/veg.spp.RDS')
veg.data.veg <-  readRDS('data/veg.data.veg.RDS')
veg.data.trans <-  readRDS('data/veg.data.trans.RDS')
# 
# veg.plot <- soilDB::get_vegplot_from_NASIS_db()
# veg.data <- soilDB::get_veg_data_from_NASIS_db()
# veg.loc <- soilDB::get_vegplot_location_from_NASIS_db()
# veg.spp <- soilDB::get_vegplot_species_from_NASIS_db()
# veg.data.veg <- veg.data$veg
# veg.data.trans <- veg.data$vegtransect
# 
# saveRDS(veg.plot, 'data/veg.plot.RDS')
# saveRDS(veg.data, 'data/veg.data.RDS')
# saveRDS(veg.loc, 'data/veg.loc.RDS')
# saveRDS(veg.spp, 'data/veg.spp.RDS')
# saveRDS(veg.data.veg, 'data/veg.data.veg.RDS')
# saveRDS(veg.data.trans, 'data/veg.data.trans.RDS')


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

colnames(veg.spp)
x <- subset(veg.spp, select = c(vegplotiid, plantsym, plantsciname, planttypegroup, plantnativity, plantheightcllowerlimit, plantheightclupperlimit, livecanopyhtbottom, livecanopyhttop, 
                                 overstorydbhmin, overstorydbhmax, speciestraceamtflag,
                                 speciescancovpct, speciescancovclass, speciescomppct, speciescompbywtpct,
                                 akstratumcoverclass, akstratumcoverclasspct)) 

unique(x$speciescancovclass)

x <- x %>% mutate(cover = case_when(
  !is.na(akstratumcoverclasspct) ~ akstratumcoverclasspct,
  !is.na(speciescancovpct) ~ speciescancovpct, + ifelse(speciestraceamtflag,0.2,0),
  speciescancovclass %in% "trace" ~ 0.05,
  speciescancovclass %in% "0.1 to 1%" ~ (0.1+1)/2,
  speciescancovclass %in% "1.1 to 2%" ~ (1+2)/2,
  speciescancovclass %in% "2 to 5%" ~ (2+5)/2,
  speciescancovclass %in% "6 to 10%" ~ (5+10)/2,
  speciescancovclass %in% "11 to 25%" ~ (10+25)/2,
  speciescancovclass %in% "26 to 50%" ~ (25+50)/2,
  speciescancovclass %in% "51 to 75" ~ (50+75)/2,
  speciescancovclass %in% "76 to 95%" ~ (75+95)/2,
  speciescancovclass %in% "> 95%" ~ (95+100)/2),
  TRUE ~ NA)
  




x <- 1:50
case_when(
  x %% 35 == 0 ~ "fizz buzz",
  x %% 5 == 0 ~ "fizz",
  x %% 7 == 0 ~ "buzz",
  TRUE ~ as.character(x)
)
                                              ))


