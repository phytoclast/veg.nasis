library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#identify location of geodatabase
mi = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb')
mu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb')
co = get_component_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb')
comonth <- sf::st_read('D:/GIS/SOIL/2021/gSSURGO_MI.gdb', 'comonth')

#extract and concatenate component table
s <-  site(mi)
h <- horizons(mi)
comonth <- comonth |> mutate(flood = ifelse(flodfreqcl %in% c("Frequent","Occasional","Rare","Very rare"), 1,0)) |> group_by(cokey) |> summarise(floodfrq=max(flood)) 


h <- h |> mutate(h50 = ifelse(hzdepb_r > 50, 50, hzdepb_r) - ifelse(hzdept_r > 50, 50, hzdept_r),
                 h150 = ifelse(hzdepb_r > 150, 150, hzdepb_r) - ifelse(hzdept_r > 150, 150, hzdept_r),
                 pH = ifelse(is.na(ph1to1h2o_r), ph01mcacl2_r, ph1to1h2o_r),
                 Bh = ifelse(grepl('Bh', hzname),1, 0),
                 rock = ifelse(grepl('Cr|R', hzname),hzdept_r, 250),
                 carb = ifelse(caco3_r > 0, hzdept_r, 250))
soil <- h |> group_by(cokey) |> summarise(sand50 = weighted.mean(sandtotal_r, h50, na.rm = T),
                                       sand150 = weighted.mean(sandtotal_r, h150, na.rm = T),
                                       pH50 = weighted.mean(pH, h50, na.rm = T),
                                       OM150 = weighted.mean(om_r, h150, na.rm = T),
                                       rockdepth = min(rock, na.rm = T),
                                       carbdepth = min(carb, na.rm = T),
                                       Bh = max(Bh, na.rm = T))
soil <- subset(mu, select=c(muname, mukey)) |> left_join(subset(s, select=c(mukey, cokey, comppct_r, majcompflag, compname, taxclname, taxorder, taxsubgrp, hydricrating, drainagecl, geomdesc))) |> left_join(soil, by=join_by(cokey==cokey)) |> left_join(comonth, by=join_by(cokey==cokey)) 


soil <- soil |> mutate(
  moist = case_when(tolower(hydricrating) %in% 'yes'~'wet',
                    drainagecl %in% c("Somewhat poorly drained","Very poorly drained","Moderately well drained","Poorly drained") ~ 'moist',
                    TRUE ~ 'dry'),
    text = case_when(rockdepth < 100 ~ 'rocky',
                   sand50 >= 70 & sand150 >= 80 | sand50 >= 80 ~ 'sandy',
                   grepl('ist', tolower(taxsubgrp)) | grepl('ist', tolower(taxorder)) |  OM150 > 20 ~ 'mucky',
                   TRUE ~ 'loamy'),
  chem = case_when(carbdepth < 50 & pH50 > 5.5 ~ 'calcareous',
                   grepl('spod',tolower(taxorder)) & (Bh > 0 | grepl('typic',tolower(taxsubgrp))) & moist %in% 'dry' & text %in% 'sandy' ~ 'spodic',
                   grepl('spod',tolower(taxorder)) & moist %in% 'dry'& text %in% 'sandy' ~ 'entic' ,
                   grepl('ult|ods',tolower(taxsubgrp)) | grepl('dys', tolower(taxclname)) ~ 'dysic',
                   grepl('alf|oll',tolower(taxsubgrp)) ~ 'euic',
                   pH50 <= 5.5 ~ 'dysic',
                   TRUE ~ 'euic'),
  flood = case_when(grepl('flood', muname) | (grepl('flood', geomdesc) & grepl('fluv', tolower(taxsubgrp))) | floodfrq > 0 ~ 'flood',
                    TRUE ~ ''))


  


#look at what the spatial layer is named
mi.layers<- sf::st_layers('D:/GIS/SOIL/2021/gSSURGO_MI.gdb')

#extract map unit polygon layer
mi.poly <- sf::st_read('D:/GIS/SOIL/2021/gSSURGO_MI.gdb', 'MUPOLYGON')

#identify map units of interest
drysand = subset(soil, text %in% 'sandy' & moist %in% 'dry' & chem %in% c('dysic','euic') & flood %in% c('') & majcompflag %in% 'Yes')
moistsand = subset(soil, text %in% 'sandy' & moist %in% 'moist' & flood %in% c('') & majcompflag %in% 'Yes')

#filter to map units of interest
drysand.poly <- subset(mi.poly, MUKEY %in% drysand$mukey)
moistsand.poly <- subset(mi.poly, MUKEY %in% moistsand$mukey)

#display map of map units
mapview(drysand.poly)

#export as shapefile
st_write(drysand.poly, 'temp/drysand.poly.shp')
st_write(moistsand.poly, 'temp/moistsand.poly.shp')

#huronstands

stands <- st_read('D:/GIS/hmnf/HNF_Stands_5_13_2022.shp')
stands <- stands |> mutate(id = SETTING_ID, veg = EV_COMMON_, year = ifelse(YEAR_OF_OR == 0, NA, YEAR_OF_OR), age2023 = 2023-year, si = SITE_INDEX, si.sp = SITE_IND_1)
jpstands <- subset(stands, grepl('jack',tolower(veg)), select=c(id, veg, year, age2023, si, si.sp))
st_write(jpstands, 'temp/jpstands.shp', append = F)
