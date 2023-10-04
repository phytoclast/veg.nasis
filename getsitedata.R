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
brick <- c(dem,slope,aspect,tpi)
site.sf <- st_as_sf(veg.site, coords	= c(x='Longitude', y='Latitude'), crs='EPSG: 4326')
site.vect <- vect(site.sf) |> project(dem)
plot(site.vect)
elevextract <- terra::extract(brick, site.vect)

site.sf.ex <- site.sf |> cbind(elevextract)
site.sf.ex <- site.sf.ex |> subset(State %in% c('Michigan','Indiana','Ohio','Illinois'), select=c(State, County,Observation_ID,  Observation_Label, elev,slope,aspect,tpi))

library(soilDB)

site.sf.soil <- soilDB::SDA_spatialQuery(site.sf)


write.csv(site.sf.ex, 'site.sf.ex.csv', row.names = F, na = "")
colnames(site.sf.ex)

