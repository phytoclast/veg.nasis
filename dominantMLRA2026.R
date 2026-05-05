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

excludedlandtypes <- c('No Digital Data Available', 'Area not surveyed, access denied', 'Water', 'Water (Inland)', 'Water (Lake Superior)', 'Water bodies', 'Water, ocean', 'Water, Saline', 'Waters of the Atlantic Ocean','Waters of the Gulf of Mexico','Waters of the Gulf of america', 'Miscellaneous water', 'Dumps, sawdust', 'Dumps, sanitary landfill','Dumps, mine', 'Dumps','Dumps, limestone', 'Sand pits','Pits, sand and gravel','Pits, quarry','Pits, quarries','Quarry','Borrow pits','Clay pit','Borrow pits and clay pits','Dumps and Pits, mine', 'Dumps, landfill','Cut and fill land','Fill land, dredged material','Fill land, sandy', 'Strip mine', 'Dumps-Pits complex, mined land', 'Gravel and sand pits','Mine pits','Pits','Pits and dumps, mine','Pits, borrow','Pits, gravel and sand','Mine pits','Pits-Dumps complex','Urban land','Mine dumps, coal','Strip mines, acid','Strip mine spoil','Sewage lagoon','Sewage lagoons','Settling pond','Dam','Large dam','Sanitary landfill','Cut and fill land','Fill land','Landfill')

#download soil geodatabase for a state https://nrcs.app.box.com/v/soils/folder/233398887779
#Name of folder where all zipped geodatabases are found (modify for your own system)
soilpath <- 'C:/GIS/SOIL/2025'
MLRA <- st_read('C:/GIS/Ecoregion/MLRA_2018.shp')

#Name of state to extract from (insert your own two letter state code or create a loop to aggregate multiple states)
soilstates <- c('CT','DC','DE','IL','IN','KY','MA','MD','ME','MI','NH','NJ','NY','OH','PA','RI','TN','VA','VT','WI','WV','NC')
#download soil legends by state
for(i in 1:length(soilstates)){#i=17
  soilstate <- soilstates[i]
  lmufilename <- paste0('gdbtemp/lmu_',soilstate,'.RDS')
  if(!file.exists(lmufilename)){
    lmu <- get_mapunit_from_NASISWebReport(areasymbol=paste0(soilstate,'%'))
    saveRDS(lmu,lmufilename)}}



#create 90-m rasters by state with water and artificial surfaces excluded
for(i in 1:length(soilstates)){#i=17
  soilstate <- soilstates[i]
  #path to the zip file of interest
  filelocation <- paste0(soilpath,'/gSSURGO_',soilstate,'.zip')
  rastfilename <- paste0('gdbtemp/mukeys90m_',soilstate,'.tif')
  lmufilename <- paste0('gdbtemp/lmu_',soilstate,'.RDS')
  #name of file at unzipped location
  dezip <- paste0('gdbtemp/gSSURGO_',soilstate,'.gdb')
  if(!file.exists(dezip)){
    unzip(filelocation, exdir ='gdbtemp')}
  #s.poly <- sf::st_read(dezip, 'MUPOLYGON')
  #use this for raster feature
  if(!file.exists(rastfilename)){
    lmu <- readRDS(lmufilename)
    excludmu <- lmu |> subset(!tolower(muname) %in% tolower(excludedlandtypes))
    excludmu <- excludmu[,c('lmapunitiid','lmapunitiid')]|> as.matrix()
    
    s.rast <- terra::rast(dezip, 'MURASTER_10m') 
    s.excl <- classify(s.rast, excludmu, others=NA)
    s.rast90 <- aggregate(s.excl,fun='modal', fact=9, filename=rastfilename, overwrite=TRUE, na.rm=TRUE)
  }}





#now either merge together and rasterize the MLRA layer
#Zonal statistics 
#
soilstate <- 'MI'
rastfilename <- paste0('gdbtemp/mukeys90m_',soilstate,'.tif')
s.rast90 <- rast(rastfilename)

r.mlra <- rasterize(vect(MLRA), s.rast90, field='LRU')
plot(r.mlra)
r <- c(s.rast90,r.mlra)
zstats <- terra::crosstab(r, long=TRUE)






writeRaster(hydric,'gdbtemp/mitemp.tif')
rvar <- subset(lmu,select=c('lmapunitiid','pct_hydric')) |> as.matrix()
hydric <- classify(s.rast90, rvar, others=NA)
plot(hydric)