setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
library(sf)
library(terra)
library(soilDB)
library(aqp)

veg.spp <- read.delim('data/Observed_Species.txt')
veg.site <- read.delim('data/Sites.txt')
mlra <- st_read("C:/a/Ecological_Sites/GIS/Ecoregion/RegionalReview_296/CONUS2.shp")
dem <- rast('D:/scripts/R12W/dem.tif'); names(dem) <- 'elev'
slope <- rast('D:/scripts/R12W/slope.tif'); names(slope) <- 'slope'
aspect <- rast('D:/scripts/R12W/aspect.tif'); names(aspect) <- 'aspect'
tpi <- rast('D:/scripts/R12W/tpi.tif'); names(tpi) <- 'tpi'
water <- read_sf('C:/a/Ecological_Sites/GIS/Base/water.shp')
veg.site <- subset(veg.site,Latitude != 0 & Observer_Code %in% c('BEL.JH', 'TOL.NB', 'GRR.NJL', 'GRR.GJS') &
                     Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #
rownames(veg.site) <- NULL
# nat500 <- rast('D:/scripts/demtools/lc.mean.tif')|> project(dem)|> extend(dem); names(nat500) <- 'nat500' 
# writeRaster(nat500, 'D:/scripts/R12W/nat500.tif', overwrite=T)
nat500 <- rast('D:/scripts/R12W/nat500.tif')
brick <- c(dem,slope,aspect,tpi,nat500)
site.sf <- st_as_sf(veg.site, coords	= c(x='Longitude', y='Latitude'), crs='EPSG: 4326')
water <- water |> subset(km2 >= 2e+05) |> vect()  |> project(dem)
site.vect <- vect(site.sf) |> project(dem)
water.dist <- distance(site.vect, water, pairwise=FALSE, unit="m")
colnames(water.dist) <- c('ocean','lake')
veg.site <- cbind(veg.site, water.dist)
plot(site.vect)
elevextract <- terra::extract(brick, site.vect)

site.sf.ex <- site.sf |> cbind(elevextract)
site.sf.ex <- site.sf.ex |> subset(select=c(State, County,Observation_ID,  Observation_Label, elev,slope,aspect,tpi, nat500)) |> st_drop_geometry()
write.csv(site.sf.ex, 'sitedata/site.sf.ex.csv', row.names = F)
#get mlra

site.mlra <- site.sf |> st_transform(crs(mlra)) |> st_intersection(mlra)
site.mlra <- subset(site.mlra, select=c(Observation_ID, LRU)) |> st_drop_geometry()
write.csv(site.mlra, 'sitedata/site.mlra.csv', row.names = F)


# library(soilDB)
# ssurgo=NULL
# ssurgo.i<-NULL
# for (i in 1:nrow(site.sf)){#i=400 length(db)
#   obs.geom.i <- site.sf[i,]
#   ssurgo.i <- SDA_spatialQuery(
#     obs.geom.i,
#     what = "mukey",
#     geomIntersection = FALSE,
#     db = c("SSURGO")
#   )
#   if(class(ssurgo.i)[1] %in% "try-error"){ssurgo.i <- cbind(mukey="try-error", muname="try-error")}
#   if(!is.null(ssurgo.i)){
#     ssurgo.i <-  cbind(list(obs.id = obs.geom.i$Observation_ID[1]), ssurgo.i) %>% as.data.frame()}
#   else{
#     ssurgo.i <-  list(obs.id = obs.geom.i$Observation_ID[1], mukey = NA, muname = NA)}
#   if(is.null(ssurgo)){
#     ssurgo <- ssurgo.i
#   }else{
#     ssurgo <- rbind(ssurgo,ssurgo.i)
#   }
# }
# saveRDS(ssurgo, 'sitedata/ssurgo.RDS')
# ssurgo <- readRDS('sitedata/ssurgo.RDS')

# #Alternative extract directly from GDB ----
# mi.poly <- sf::st_read('D:/GIS/SOIL/2021/gSSURGO_MI.gdb', 'MUPOLYGON')
# oh.poly <- sf::st_read('D:/GIS/SOIL/2021/gSSURGO_OH.gdb', 'MUPOLYGON')
# in.poly <- sf::st_read('D:/GIS/SOIL/2021/gSSURGO_IN.gdb', 'MUPOLYGON')
# il.poly <- sf::st_read('D:/GIS/SOIL/2021/gSSURGO_IL.gdb', 'MUPOLYGON')
# nj.poly <- sf::st_read('D:/GIS/SOIL/2021/gSSURGO_NJ.gdb', 'MUPOLYGON')
# all.poly <- rbind(mi.poly,oh.poly,in.poly,il.poly,nj.poly)
# 
# site.sfmu <- site.sf |> st_transform(crs(all.poly))
# 
# site.inter <- sf::st_intersection(site.sfmu, all.poly)
# colnames(site.inter)
# site.dropgeo <- site.inter |> st_drop_geometry() |> subset(select = c(Observation_ID, MUKEY))
# saveRDS(site.dropgeo, 'sitedata/site.dropgeo.RDS')
site.dropgeo <- readRDS('sitedata/site.dropgeo.RDS')
#Soil component tables----
library(soilDB)
library(aqp)
# remotes::install_github('ncss-tech/soilDB')
#
# source('fetchGDB.R')
#
# MI = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb', childs = TRUE)
# IN = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IN.gdb', childs = TRUE)
# OH = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_OH.gdb', childs = TRUE)
# IL = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IL.gdb', childs = TRUE)
# NJ = fetchGDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_NJ.gdb', childs = TRUE)
# 
# MImu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_MI.gdb')
# INmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IN.gdb')
# OHmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_OH.gdb')
# ILmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_IL.gdb')
# NJmu = get_mapunit_from_GDB(dsn = 'D:/GIS/SOIL/2021/gSSURGO_NJ.gdb')
# 
# ssurgotabs <- c(MI,IN,OH,IL,NJ) |> subset(majcompflag %in% 'Yes')
# ssurgomu <- rbind(MImu,INmu,OHmu,ILmu,NJmu)
# saveRDS(ssurgotabs, 'ssurgo/ssurgotabs.RDS')
# saveRDS(ssurgomu, 'ssurgo/ssurgomu.RDS')
ssurgotabs <- readRDS('ssurgo/ssurgotabs.RDS')
ssurgomu <- readRDS('ssurgo/ssurgomu.RDS')

sssites <- site(ssurgotabs) |> left_join(ssurgomu)
sshorz <- horizons(ssurgotabs)

#Pedon data ---- 
veg.site.mi <- subset(veg.site, State %in% c('Michigan', 'Illinois', 'Indiana', 'Ohio'))
max(veg.site.mi$Latitude)
min(veg.site.mi$Latitude)
max(veg.site.mi$Longitude)
min(veg.site.mi$Longitude)
veg.site.nj <- subset(veg.site, State %in% c('New Jersey'))
max(veg.site.nj$Latitude)
min(veg.site.nj$Latitude)
max(veg.site.nj$Longitude)
min(veg.site.nj$Longitude)

# nasispedons <- fetchNASIS(from = 'pedons', SS=FALSE)
# saveRDS(nasispedons, 'sitedata/nasispedons.RDS')
nasispedons <- readRDS('sitedata/nasispedons.RDS')

nassites <- site(nasispedons)
nashorz <- horizons(nasispedons)

n <- nrow(veg.site)
for (i in 1:n){#i=484
  this <- veg.site[i,]
  thislat <- this$Latitude
  thislon <- this$Longitude
  thisrecord <-  subset(this, select = c(Observation_ID, Observation_Label, User_Pedon_ID, Latitude, Longitude, State, County, ocean, lake))
  nassites <- nassites |> mutate(distance = (((thislat - latstddecimaldegrees)/360*40041.47*1000)^2 +
                             ((thislon - longstddecimaldegrees)/360*40041.47*1000*cos(thislat/2/360*2*3.141592))^2)^0.5)
  mindist <- min(nassites$distance, na.rm = T)
  thispedon <- subset(nassites, distance <= mindist, select=c(pedon_id, distance))[1,]
  if(i==1){myrecords <- cbind(thisrecord, thispedon)}else{
    myrecords <- rbind(myrecords, cbind(thisrecord, thispedon))
  }}

nasss <- subset(nashorz, select=c(texture, sand, clay))
nas.textures <- texcl_to_ssc(texcl=nashorz$texture)
colnames(nas.textures) <- c("sand1", "silt1", "clay1")
nashorz <- nashorz |> cbind(nas.textures)
nashorz <- nashorz |> mutate(sand = ifelse(is.na(sand), sand1, sand))

#back to processing the ssurgo data ----
sshorzsand <- sshorz |> mutate(th050 = ifelse(hzdepb_r > 50,50,hzdepb_r)-ifelse(hzdept_r > 50,50,hzdept_r),
                               th150 = ifelse(hzdepb_r > 150,150,hzdepb_r)-ifelse(hzdept_r > 150,150,hzdept_r)) |>
  group_by(cokey) |> summarise(sand050 = sum(sandtotal_r*th050, na.rm = TRUE)/sum(th050, na.rm = TRUE),
                               sand150 = sum(sandtotal_r*th150, na.rm = TRUE)/sum(th150, na.rm = TRUE))
sshorzpH <- sshorz |> mutate(th050 = ifelse(hzdepb_r > 50,50,hzdepb_r)-ifelse(hzdept_r > 50,50,hzdept_r),
                             pH = ifelse(!is.na(ph1to1h2o_r), ph1to1h2o_r, ph01mcacl2_r)) |> subset(!is.na(pH)) |>
  group_by(cokey) |> summarise(pH50 = sum(pH*th050, na.rm = TRUE)/sum(th050, na.rm = TRUE))
#back to processing the pedon data ----
nashorzsand <- nashorz |> mutate(th050 = ifelse(hzdept > 50,50,hzdepb)-ifelse(hzdept > 50,50,hzdept),
                               th150 = ifelse(hzdepb > 150,150,hzdepb)-ifelse(hzdept > 150,150,hzdept)) |>
  group_by(peiid) |> summarise(sand050 = sum(sand*th050, na.rm = TRUE)/sum(th050, na.rm = TRUE),
                               sand150 = sum(sand*th150, na.rm = TRUE)/sum(th150, na.rm = TRUE))
nashorzpH <- nashorz |> mutate(th050 = ifelse(hzdept > 50,50,hzdepb)-ifelse(hzdept > 50,50,hzdept),
                             pH = ifelse(!is.na(phfield), phfield, ifelse(!is.na(effclass) & !effclass %in% 'none',8,NA))) |> subset(!is.na(pH)) |>
  group_by(peiid) |> summarise(pH50 = sum(pH*th050, na.rm = TRUE)/sum(th050, na.rm = TRUE))

#back to processing the ssurgo data ----
sshorzrock <- sshorz |> mutate(rock = ifelse(texture %in% c('UWB','WB','BR')|grepl('Cr', hzname), 'BR','Soil')) |> subset(rock %in% 'BR') |> group_by(cokey) |> summarise(rockdepth = min(hzdepb_r))

sssitesflood <- sssites |> mutate(flood = ifelse(muname %in% c("flooded") | compname %in% c("Alluvial land")|
                                                  (grepl("flood",tolower(geomdesc)) & (grepl("fluv",tolower(taxsubgrp))|grepl("psam",tolower(taxsubgrp))|grepl("cumu",tolower(taxsubgrp)))),1,0))

#back to processing the pedon data ----
nassitesplus <- nassites |> subset(select= c(obs_date, siteiid, peiid, pedon_id, ecositeid, drainagecl, taxonname, taxclname, taxsubgrp, bedrckdepth, flodfreqcl, landform_string)) |> left_join(nashorzsand) |> left_join(nashorzpH)


#back to processing the ssurgo data ----
ssitesplus <- sssites |> subset(select= c(mukey, cokey, drainagecl,comppct_r, compname, taxclname, taxsubgrp)) |> left_join(sshorzsand)|> left_join(sshorzrock)|> left_join(sssitesflood) |> left_join(sshorzpH) |> mutate(mukey = as.integer(mukey))

# ssurgoplus <- ssurgo |> left_join(ssitesplus) # previously used link to online ssurgo
ssurgoplus <- site.dropgeo |> mutate(mukey = as.numeric(MUKEY), MUKEY = NULL) |> left_join(ssitesplus, by=join_by(mukey==mukey)) #join to mukey from downloaded gdb 

ssurgoplus <- ssurgoplus |> left_join(veg.site[,c('Observation_ID', 'ocean','lake')], by=join_by(Observation_ID == Observation_ID))

myrecordsplus <- myrecords |> left_join(nassitesplus, by=join_by(pedon_id==pedon_id))

#assess flood,sand,bedrock,wetland statuses of site, then of current condition
unique(ssurgoplus$drainagecl)
unique(myrecordsplus$drainagecl)
ssurgoplus <- ssurgoplus |> 
  mutate(wet = ifelse(hydricrating %in% "Yes" | drainagecl %in% c("Poorly drained","Very poorly drained"),1,0),
                                   moist = ifelse(drainagecl %in% c("Moderately well drained","Somewhat poorly drained"),1,0),
                                   dry = ifelse(wet+moist > 0,0,1),
                                   rockdepth = ifelse(is.na(rockdepth), 500,rockdepth), 
                                   mucky = ifelse(grepl('istic', taxsubgrp) | taxorder %in% 'Histosols', 1,0),
                                   sandy = ifelse(sand050 >= 80 | (sand050 >= 70 & sand150 >= 80), 1,0),
                                   rock = ifelse((grepl('Lithic', taxsubgrp))|rockdepth <= 100,1,0),
                                   coastal = ifelse(ocean <= 1000| lake <= 500, 1,0),
                                   salty = ifelse(grepl('tidal',landform), 1,0),
         dysic = ifelse(is.na(pH50), 
                        ifelse(grepl('dys',taxclname)| grepl('ult',taxclname)| grepl('ods',taxclname),1,0),
                        ifelse(pH50 <=5.5,1,0)),
         euic = ifelse(is.na(pH50), 
                        ifelse(grepl('eu',taxclname)| grepl('alf',taxclname)| grepl('oll',taxclname),1,0),
                        ifelse(pH50 >5.5,1,0)))

myrecordsplus <- myrecordsplus |> mutate(p_wet = ifelse(drainagecl %in% c("poorly","very poorly"),1,0),
                                   p_moist = ifelse(drainagecl %in% c("moderately well","somewhat poorly"),1,0),
                                   p_dry = ifelse(p_wet+p_moist > 0,0,1),
                                   p_mucky = ifelse(grepl('ist', taxsubgrp), 1,0),
                                   bedrckdepth=ifelse(is.na(bedrckdepth),500,bedrckdepth),
                                   p_sandy = ifelse(sand050 >= 80 | (sand050 >= 70 & sand150 >= 80), 1,0),
                                   p_rock = ifelse((grepl('Lithic', taxsubgrp))|bedrckdepth <= 100,1,0),
                                   p_coastal = ifelse(ocean <= 1000| lake <= 500, 1,0),
                                   p_flood = ifelse((!flodfreqcl %in% 'none' & !is.na(flodfreqcl)) | grepl('flood',landform_string),1,0),
                                   p_salty = ifelse(pedon_id %in% c('2016MI037001','2016MI037002') | grepl('tidal',landform_string), 1,0),
                                   p_dysic = ifelse(is.na(pH50), 
                                                  ifelse(grepl('dys',taxclname)| grepl('ult',taxclname)| grepl('ods',taxclname),1,0),
                                                  ifelse(pH50 <=5.5,1,0)),
                                   p_euic = ifelse(is.na(pH50), 
                                                 ifelse(grepl('eu',taxclname)| grepl('alf',taxclname)| grepl('oll',taxclname),1,0),
                                                 ifelse(pH50 >5.5,1,0))
                                   )

write.csv(myrecordsplus, 'sitedata/myrecordsplus.csv', row.names = F)
write.csv(ssurgoplus, 'sitedata/ssurgoplus.csv', row.names = F)


#process site data  ----

myrecordsplus <- read.csv('sitedata/myrecordsplus.csv')
ssurgoplus <- read.csv('sitedata/ssurgoplus.csv')
site.sf.ex <- read.csv('sitedata/site.sf.ex.csv')
site.mlra <- read.csv('sitedata/site.mlra.csv')


colnames(ssurgoplus)

ssurgoplus1 <- ssurgoplus |> subset(select = c("Observation_ID", "mukey", "muname", "flood","wet","moist","dry" ,"mucky","sandy","rock","coastal","salty","dysic","euic"))
ssurgoplus1 <- ssurgoplus1 |> group_by(Observation_ID, mukey, muname) |> summarise(across(c("flood","wet","moist","dry" ,"mucky","sandy","rock","coastal","salty","dysic","euic"), mean)) 

myrecordsplus1 <- myrecordsplus |> subset(distance < 100, select=c("Observation_ID", "Observation_Label", "pedon_id","distance","taxonname","taxclname","p_flood","p_wet","p_moist","p_dry" ,"p_mucky","p_sandy","p_rock","p_coastal","p_salty","p_dysic","p_euic")) 

ssurgoplus1 <- ssurgoplus1 |> left_join(myrecordsplus1)
ssurgoplus1 <- ssurgoplus1 |> mutate(
flood = ifelse(is.na(p_flood),flood,p_flood), p_flood=NULL,
 wet = ifelse(is.na(p_wet),wet,p_wet), p_wet=NULL,
 moist = ifelse(is.na(p_moist),moist,p_moist), p_moist=NULL,
 dry = ifelse(is.na(p_dry),dry,p_dry), p_dry=NULL,
 mucky = ifelse(is.na(p_mucky),mucky,p_mucky), p_mucky=NULL,
 sandy = ifelse(is.na(p_sandy),sandy,p_sandy), p_sandy=NULL,
 rock = ifelse(is.na(p_rock),rock,p_rock), p_rock=NULL,
 coastal = ifelse(is.na(p_coastal),coastal,p_coastal), p_coastal=NULL,
 salty = ifelse(is.na(p_salty),salty,p_salty), p_salty=NULL,
 dysic = ifelse(is.na(p_dysic),dysic,p_dysic), p_dysic=NULL,
 euic = ifelse(is.na(p_euic),euic,p_euic), p_euic=NULL
)
