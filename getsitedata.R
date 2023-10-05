setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
library(sf)
library(terra)
veg.spp <- read.delim('data/Observed_Species.txt')
veg.site <- read.delim('data/Sites.txt')
dem <- rast('D:/scripts/R12W/dem.tif'); names(dem) <- 'elev'
slope <- rast('D:/scripts/R12W/slope.tif'); names(slope) <- 'slope'
aspect <- rast('D:/scripts/R12W/aspect.tif'); names(aspect) <- 'aspect'
tpi <- rast('D:/scripts/R12W/tpi.tif'); names(tpi) <- 'tpi'

veg.site <- subset(veg.site,Latitude != 0 & Observer_Code %in% c('BEL.JH', 'TOL.NB', 'GRR.NJL', 'GRR.GJS') & 
                     Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #


brick <- c(dem,slope,aspect,tpi)
site.sf <- st_as_sf(veg.site, coords	= c(x='Longitude', y='Latitude'), crs='EPSG: 4326')
site.vect <- vect(site.sf) |> project(dem)
plot(site.vect)
elevextract <- terra::extract(brick, site.vect)

site.sf.ex <- site.sf |> cbind(elevextract)
site.sf.ex <- site.sf.ex |> subset(select=c(State, County,Observation_ID,  Observation_Label, elev,slope,aspect,tpi))

library(soilDB)
ssurgo=NULL
ssurgo.i<-NULL
for (i in 1:nrow(site.sf)){#i=400 length(db)
  obs.geom.i <- site.sf[i,]
  ssurgo.i <- SDA_spatialQuery(
    obs.geom.i,
    what = "mukey",
    geomIntersection = FALSE,
    db = c("SSURGO")
  )
  if(class(ssurgo.i)[1] %in% "try-error"){ssurgo.i <- cbind(mukey="try-error", muname="try-error")}
  if(!is.null(ssurgo.i)){
    ssurgo.i <-  cbind(list(obs.id = obs.geom.i$Observation_ID[1]), ssurgo.i) %>% as.data.frame()}
  else{
    ssurgo.i <-  list(obs.id = obs.geom.i$Observation_ID[1], mukey = NA, muname = NA)}
  if(is.null(ssurgo)){
    ssurgo <- ssurgo.i
  }else{
    ssurgo <- rbind(ssurgo,ssurgo.i)
  }
}
saveRDS(ssurgo, 'sitedata/ssurgo.RDS')

library(soilDB)
library(aqp)
MI = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb', childs = TRUE)
IN = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IN.gdb', childs = TRUE)
OH = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_OH.gdb', childs = TRUE)
IL = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IL.gdb', childs = TRUE)
NJ = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_NJ.gdb', childs = TRUE)

MImu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb')
INmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IN.gdb')
OHmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_OH.gdb')
ILmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IL.gdb')
NJmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_NJ.gdb')

# MIco = get_component_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb', childs = TRUE)
# INco = get_component_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IN.gdb')
# OHco = get_component_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_OH.gdb')
# ILco = get_component_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IL.gdb')
# NJco = get_component_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_NJ.gdb')


  

ssurgotabs <- c(MI,IN,OH,IL,NJ) |> subset(majcompflag %in% 'Yes')
ssurgomu <- rbind(MImu,INmu,OHmu,ILmu,NJmu)
# ssurgoco <- rbind(MIco,INco,OHco,ILco,NJco)
sssites <- site(ssurgotabs) |> left_join(ssurgomu)
sshorz <- horizons(ssurgotabs)

# veg.site.mi <- subset(veg.site, State %in% c('Michigan', 'Illinois', 'Indiana', 'Ohio'))
# max(veg.site.mi$Latitude)
# min(veg.site.mi$Latitude)
# max(veg.site.mi$Longitude)
# min(veg.site.mi$Longitude)
# veg.site.nj <- subset(veg.site, State %in% c('New Jersey'))
# max(veg.site.nj$Latitude)
# min(veg.site.nj$Latitude)
# max(veg.site.nj$Longitude)
# min(veg.site.nj$Longitude)
# 
# nasispedons <- fetchNASIS(from = 'pedons', SS=FALSE)
# saveRDS(nasispedons, 'sitedata/nasispedons.RDS')
nasispedons <- readRDS('sitedata/nasispedons.RDS')

nassites <- site(nasispedons)
nashorz <- horizons(nasispedons)

n <- nrow(veg.site)
for (i in 1:n){#i=1
  this <- veg.site[i,]
  thislat <- this$Latitude
  thislon <- this$Longitude
  thisrecord <-  subset(this, select = c(Observation_ID, Observation_Label, User_Pedon_ID, Latitude, Longitude, State, County))
  nassites <- nassites |> mutate(distance = (((thislat - latstddecimaldegrees)/360*40041.47*1000)^2 +
                             ((thislon - longstddecimaldegrees)/360*40041.47*1000*cos(thislat/2/360*2*3.141592))^2)^0.5)
  mindist <- min(nassites$distance, na.rm = T)
  thispedon <- subset(nassites, distance <= mindist, select=c(pedon_id, distance))[1,]
  if(i==1){myrecords <- cbind(thisrecord, thispedon)}else{
    myrecords <- rbind(myrecords, cbind(thisrecord, thispedon))
  }}
  


sshorzsand <- sshorz |> mutate(th050 = ifelse(hzdepb_r > 50,50,hzdepb_r)-ifelse(hzdept_r > 50,50,hzdept_r),
                               th150 = ifelse(hzdepb_r > 150,150,hzdepb_r)-ifelse(hzdept_r > 150,150,hzdept_r)) |> 
  group_by(cokey) |> summarise(sand050 = sum(sandtotal_r*th050, na.rm = TRUE)/sum(th050, na.rm = TRUE),
                               sand150 = sum(sandtotal_r*th150, na.rm = TRUE)/sum(th150, na.rm = TRUE))

nashorzsand <- nashorz |> mutate(th050 = ifelse(hzdept > 50,50,hzdepb)-ifelse(hzdept > 50,50,hzdept),
                               th150 = ifelse(hzdepb > 150,150,hzdepb)-ifelse(hzdept > 150,150,hzdept)) |> 
  group_by(peiid) |> summarise(sand050 = sum(sand*th050, na.rm = TRUE)/sum(th050, na.rm = TRUE),
                               sand150 = sum(sand*th150, na.rm = TRUE)/sum(th150, na.rm = TRUE))

sshorzrock <- sshorz |> mutate(rock = ifelse(texture %in% c('UWB','WB','BR')|grepl('Cr', hzname), 'BR','Soil')) |> subset(rock %in% 'BR') |> group_by(cokey) |> summarise(rockdepth = min(hzdepb_r))

# sssitesflood <- sssites |> mutate(flood = ifelse(muname %in% c("flooded") | compname %in% c("Alluvial land")|
#                                                    (grepl("flood",landform_string) & (grepl("fluv",taxsubgrp)|grepl("psam",taxsubgrp)|grepl("cumu",s$taxsubgrp))),'flood','not'))

sssitesflood <- sssites |> mutate(flood = ifelse(muname %in% c("flooded") | compname %in% c("Alluvial land"),'flood','not')) |> subset(select=c(cokey, flood))


nassitesplus <- nassites |> subset(select= c(obs_date, siteiid, peiid, pedon_id, ecositeid, drainagecl, taxonname, taxclname, taxsubgrp, bedrckdepth, flodfreqcl)) |> left_join(nashorzsand)



ssitesplus <- sssites |> subset(select= c(mukey, cokey, drainagecl,comppct_r, compname, taxclname, taxsubgrp)) |> left_join(sshorzsand)|> left_join(sshorzrock)|> left_join(sssitesflood) |> mutate(mukey = as.integer(mukey))

ssurgoplus <- ssurgo |> left_join(ssitesplus)

myrecordsplus <- myrecords |> left_join(nassitesplus, by=join_by(pedon_id==pedon_id))

#assess flood,sand,bedrock,wetland statuses of site, then of current condition