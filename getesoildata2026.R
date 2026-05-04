#remotes::install_github("ncss-tech/soilDB@fix-466")
library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
library(future)
library(rgbif)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
timeA = Sys.time()
#download soil geodatabase for a state https://nrcs.app.box.com/v/soils/folder/233398887779

#Name of folder where all zipped geodatabases are found (modify for your own system)
soilpath <- 'C:/GIS/SOIL/2025'
#Name of state to extract from (insert your own two letter state code or create a loop to aggregate multiple states)
soilstate <- 'MI'

#path to the zip file of interest
filelocation <- paste0(soilpath,'/gSSURGO_',soilstate,'.zip')
unzip(filelocation, exdir ='gdbtemp')
#name of file at unzipped location
dezip <- paste0('gdbtemp/gSSURGO_',soilstate,'.gdb')
#Extract soil profile collection which aggregates several tables
mysoil = fetchGDB(dsn = dezip)
#Extract additional tables as needed. To see list of table names use:# sf::st_layers(dezip)
#mu = get_mapunit_from_GDB(dsn = dezip)#broken
#co = get_component_from_GDB(dsn = dezip)
lu <- sf::st_read(dezip, 'legend')
mu <- sf::st_read(dezip, 'mapunit')
co =  sf::st_read(dezip, 'component')
comonth <- sf::st_read(dezip, 'comonth')
coecoclass <- sf::st_read(dezip, 'coecoclass')
#Extract GIS layers
#use this for polygon feature
#s.poly <- sf::st_read(dezip, 'MUPOLYGON')
#use this for raster feature
s.rast <- terra::rast(dezip, 'MURASTER_10m') 
plot(s.rast)
#extract and concatenate component table
s <-  site(mysoil)
h <- horizons(mysoil)
comonth <- comonth |> mutate(flood = ifelse(flodfreqcl %in% c("Frequent","Occasional","Rare","Very rare"), 1,0)) |> group_by(cokey) |> summarise(floodfrq=max(flood))

#mean properties weighted by horizon thickness
h <- h |> mutate(h50 = ifelse(hzdepb_r > 50, 50, hzdepb_r) - ifelse(hzdept_r > 50, 50, hzdept_r),
                 h150 = ifelse(hzdepb_r > 150, 150, hzdepb_r) - ifelse(hzdept_r > 150, 150, hzdept_r),
                 pH = ifelse(is.na(ph1to1h2o_r), ph01mcacl2_r, ph1to1h2o_r),
                 rock = ifelse(grepl('Cr|R', hzname),hzdept_r, 250),
                 carb = ifelse(!is.na(caco3_r) & caco3_r > 0, hzdept_r, 250),
                 sand = ifelse(is.na(sandtotal_r) & om_r >= 35,0,sandtotal_r),
                 silt = ifelse(is.na(silttotal_r) & om_r >= 35,0,silttotal_r),
                 clay = ifelse(is.na(claytotal_r) & om_r >= 35,0,claytotal_r))
soil <- h |> group_by(cokey) |> summarise(sand50 = weighted.mean(sand, h50, na.rm = T),
                                          sand150 = weighted.mean(sand, h150, na.rm = T),
                                          clay150 = weighted.mean(clay, h150, na.rm = T),
                                          pH50 = weighted.mean(pH, h50, na.rm = T),
                                          OM150 = weighted.mean(om_r, h150, na.rm = T),
                                          rockdepth = min(rock, na.rm = T),
                                          carbdepth = min(carb, na.rm = T))
soil <- subset(mu, select=c(muname, mukey)) |> left_join(subset(s, select=c(mukey, cokey, comppct_r, compname, majcompflag, taxclname, taxorder, taxsubgrp, hydricrating, drainagecl, geomdesc))) |> left_join(soil, by=join_by(cokey==cokey)) |> left_join(comonth, by=join_by(cokey==cokey))
#create soil categories
soil <- soil |> mutate(
  moist = case_when(hydricrating %in% 'yes'~'wet',
                    drainagecl %in% c("Somewhat poorly drained","Very poorly drained","Moderately well drained","Poorly drained") ~ 'moist',
                    TRUE ~ 'dry'),
  rockdepth = ifelse(rockdepth > 50 & grepl('lithic',tolower(taxsubgrp)), 25, rockdepth),
  
  chem = case_when(carbdepth < 100 ~ 'calcareous',
                   grepl('ult|ods',tolower(taxsubgrp)) | grepl('dys', tolower(taxclname)) ~ 'dysic',
                   grepl('alf|oll',tolower(taxsubgrp)) ~ 'euic',
                   pH50 < 5.5 ~ 'dysic',
                   TRUE ~ 'euic'),
  text = case_when(rockdepth < 100 ~ 'rocky',
                   sand50 >= 70 & sand150 >= 80 | sand50 >= 80 ~ 'sandy',
                   grepl('ist', tolower(taxsubgrp)) | grepl('ist', tolower(taxorder)) |  OM150 > 20 ~ 'mucky',
                   TRUE ~ 'loamy'),
  flood = case_when(grepl('flood', muname) | (grepl('flood', geomdesc) & grepl('fluv', tolower(taxsubgrp))) | floodfrq > 0 ~ 'flood',
                    TRUE ~ ''),
  watertable = case_when(drainagecl %in% c('Well drained', 'Somewhat excessively drained', 'Excessively drained') ~ 250,
                         drainagecl %in% c('Moderately well drained') ~ 100,
                         drainagecl %in% c('Somewhat poorly drained') ~ 50,
                         drainagecl %in% c('Poorly drained') ~ 25,
                         drainagecl %in% c('Very poorly drained') ~ 0,
                         TRUE ~ NA),
  hydric = ifelse(hydricrating %in% 'Yes',1,0))

#summarize mapunits by major component weighted by composition
soilsummary <- soil |> subset(majcompflag %in% 'Yes' & !is.na(sand50)) |>
  group_by(mukey) |>
  summarise(sand50 = weighted.mean(sand50, comppct_r, na.rm=T),
            sand150 = weighted.mean(sand150, comppct_r, na.rm=T),
            clay150 = weighted.mean(clay150, comppct_r, na.rm=T),
            pH50 = weighted.mean(pH50, comppct_r, na.rm=T),
            OM150 = weighted.mean(OM150, comppct_r, na.rm=T),
            hydric = weighted.mean(hydric, comppct_r, na.rm=T),
            watertable = weighted.mean(watertable, comppct_r, na.rm=T),
            carbdepth = weighted.mean(carbdepth, comppct_r, na.rm=T),
            rockdepth = weighted.mean(rockdepth, comppct_r, na.rm=T),
            floodfrq = weighted.mean(floodfrq, comppct_r, na.rm=T)) |> ungroup() |>
  subset(select=c(mukey, sand50, sand150, clay150, pH50, OM150, watertable, hydric, rockdepth, carbdepth, floodfrq))

timeA = Sys.time()
s.rast90 <- aggregate(s.rast,fun='modal', fact=9)
Sys.time() - timeA 

valclass <- soilsummary[,c('mukey','sand150')] |> as.matrix()
r <- classify(s.rast90a, valclass, others=NA)
writeRaster(r, 'gdbtemp/mi-a.tif', overwrite=T)
















#GBIF
#geographic constrants of gbif lookups
lat <- paste0(as.character(40), ",",as.character(50))
lon <- paste0(as.character(-90), ",",as.character(-80))

#get example plant points to compare
taxa <- c('Carex atlantica','Carex interior')
plt <- NULL
for(i in 1:length(taxa)){
  taxon = taxa[i]
  try({   #suppress server error and move on to next taxon
    p0 <- occ_search(limit=1000,
                     phylumKey = 7707728, scientificName = taxon,
                     #gadmGid="USA.23_1", #geographic admin unit
                     decimalLatitude=lat, decimalLongitude =lon)
  })
  p0 <- p0$data 
  p0 <- p0 |> select("key","scientificName",
                     "decimalLatitude","decimalLongitude",
                     "continent","stateProvince","year","month","day","basisOfRecord",
                     "occurrenceStatus","taxonKey","name","institutionCode","recordedBy","identifiedBy")
  if(is.null(plt)){plt <- p0}else{plt <- rbind(plt,p0)}
};rm(p0)
#Get species name only
plt$binomial <- vegnasis::extractTaxon(plt$scientificName, report='binomial')
plt$lat <- plt$decimalLatitude
plt$lon <- plt$decimalLongitude
#Create point spatial features
plt.geo <- st_as_sf(plt, coords = c(x='lon', y='lat'), crs='epsg:4326')








plot(r)
#extract soil mapunits to points 

soilpts <- terra::extract(s.rast, vect(plt.geo))
plt.join <- cbind(plt, soilpts) |> mutate(mukey = as.character(MUKEY)) #|> subset(!is.na(mukey))
plt.join <- plt.join |> left_join(soilsummary, by=join_by(mukey==mukey))  |> subset(!is.na(sand50)) 


library(ggplot2)

ggplot(plt.join)+
  geom_boxplot(aes(x=binomial, y=hydric))

library(climatools)
vars = c('lat', 'lon', 'sand50', 'sand150', 'clay150', 'pH50', 'OM150', 'watertable', 'hydric', 'rockdepth', 'carbdepth', 'floodfrq')
scompare <- plt.join |> mutate(sp1 = ifelse(binomial %in% 'Carex atlantica',1,0), sp2 = ifelse(binomial %in% 'Carex interior',1,0)) |> subset(sp1+sp2 > 0)
varimportance <- climatools::find.multithreshold(scompare, class1 = 'sp1', class2 = 'sp2', variables = vars)

