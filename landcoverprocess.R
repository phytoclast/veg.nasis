library(terra)
library(dplyr)
library(sf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

focalmean <- function(x, r){
  #establish aggregating factor when radius is too large
  fc = floor(r/res(x)[1]/9+1)
  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'mean',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  fm <- focalMat(x1, d=r, type = 'circle')
  #exclude outer portion of circle and ensure max/min values are only multiplied by 1
  fm.na <- ifelse(fm > 0, 1, NA)
  x1.mean <- focal(x1, fm.na, fun='mean', na.rm=T)
  #restore resolution in result
  if(fc > 1){
    x1.mean <- project(x1.mean, x)
  }
  return(x1.mean)} 

veg.site <- read.delim('data/Sites.txt')
veg.site <- subset(veg.site,Latitude != 0 & Observer_Code %in% c('BEL.JH', 'TOL.NB', 'GRR.NJL', 'GRR.GJS') &
                     Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #
rownames(veg.site) <- NULL

site.sf <- st_as_sf(veg.site, coords	= c(x='Longitude', y='Latitude'), crs='EPSG: 4326')
veg.site.mi <- subset(site.sf, State %in% c('Michigan', 'Illinois', 'Indiana', 'Ohio'))
veg.site.nj <- subset(site.sf, State %in% c('New Jersey'))

lc <- rast('C:/a/Ecological_Sites/GIS/Vegetation/National_Landcover_l48_eslf_v3_5_2018/l48_IVC_Existing_v3_5_2018.tif')
veg.site.mi <- veg.site.mi |> st_transform(crs(lc))
veg.site.nj <- veg.site.nj |> st_transform(crs(lc))

ex1 <- ext(veg.site.nj)
ex2 <- ex1
ex2[1] <- ex2[1] -10000
ex2[2] <- ex2[2] +10000
ex2[3] <- ex2[3] -10000
ex2[4] <- ex2[4] +10000
lc.crop <-  crop(lc, ex2)
rat <- levels(lc.crop) |> as.data.frame()

lc.natural <- ifel(lc.crop %in% c(-99, 11)| is.na(lc.crop), NA, ifel(lc.crop < 20 | lc.crop > 2000,1,0))
plot(lc.natural)

writeRaster(lc.natural, 'lc.naturalnj.tif', overwrite=T)
lc.mean <- focalmean(lc.natural, 500)
writeRaster(lc.mean, 'lc.meannj.tif', overwrite=T)

