#remotes::install_github("ncss-tech/soilDB@fix-466") #This is a fix to SoilDB for fetch GDB functions
library(soilDB)
library(aqp)
library(sf)
# library(vegnasis)
library(terra)
# library(rgbif)
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
for(i in 1:length(soilstates)){#i=15
  soilstate <- soilstates[i]
  lmufilename <- paste0('gdbtemp/lmu_',soilstate,'.RDS')
  if(!file.exists(lmufilename)){
    lmu <- get_mapunit_from_NASISWebReport(areasymbol=paste0(soilstate,'%'))
    saveRDS(lmu,lmufilename)}}
# lmu1 <- get_mapunit_from_NASISWebReport(areasymbol='PA%6%')
#missing areas c(PA605, PA611, PA607, PA609, PA610)

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
    # # make borders NA to prevent interference with merging
    # excludmu <- lmu |> subset(!tolower(muname) %in% tolower(excludedlandtypes))
    # excludmu <- excludmu[,c('lmapunitiid','lmapunitiid')]|> as.matrix()
    # 
    # s.rast <- terra::rast(dezip, 'MURASTER_10m')
    # s.excl <- classify(s.rast, excludmu, others=NA)excludmu <- lmu |> subset(tolower(muname) %in% tolower(excludedlandtypes)) |> mutate(val=NA)
    #This excludes non-natural soils and water bodies so that the resample doesn't overemphasize areas that is irrelevent to MLRA identity
    excludmu <- lmu |> subset(tolower(muname) %in% tolower(excludedlandtypes)) |> mutate(val=NA)
    excludmu <- excludmu[,c('lmapunitiid','val')]|> as.matrix()
    #Get the raster from the geodatabase
    s.rast <- terra::rast(dezip, 'MURASTER_10m') 
    s.excl <- classify(s.rast, excludmu)
    #Resample to 90m
    s.rast90 <- aggregate(s.excl,fun='modal', fact=9, filename=rastfilename, overwrite=TRUE, na.rm=TRUE)
    #This removes the zeros in the border areas so that the merge can happen correctly
    s.rast90a <- ifel(s.rast90 == 0,NA,s.rast90)
    writeRaster(s.rast90a, filename=rastfilename, overwrite=TRUE)
  }}

regionalfile <- 'gdbtemp/mukeys90_NE_Region.tif'
if(!file.exists(regionalfile)){
  filez <- paste0('gdbtemp/mukeys90m_',soilstates,'.tif')
  essembledtifs <- sprc(filez)
  essembledtifs1 <- merge(essembledtifs, method="near")
  fillcracks <- focal(essembledtifs1, fun='modal', na.policy="only", na.rm=TRUE)
  writeRaster(fillcracks,regionalfile, overwrite=T)}
mkeys <- rast(regionalfile)
names(mkeys) <- 'mkeys'

r.mlra <- rasterize(vect(MLRA), mkeys, field='LRU')
plot(r.mlra)
#now either merge together and rasterize the MLRA layer
#crosstabulate statistics 
r <- c(mkeys,r.mlra)
zstats <- terra::crosstab(r, long=TRUE)

#gatherup all the mudata saved earlier
lmu <- NULL
for(i in 1:length(soilstates)){#i=17
  soilstate <- soilstates[i]
  lmufilename <- paste0('gdbtemp/lmu_',soilstate,'.RDS')
  lmu0 <- readRDS(lmufilename)
  if(is.null(lmu)){lmu=lmu0}else{lmu<-rbind(lmu,lmu0)}
  lmu0 <- NULL
}

#Crosstabulation of the MLRA overlap with mapunits
lmustat <- zstats |> left_join(lmu, by=join_by(mkeys==lmapunitiid))
#merge LRUs that are not distinguished
lmustat <- lmustat |> mutate(MLRA = case_when(grepl('139', LRU) ~ '139',
                                              grepl('111B', LRU) ~ '111B',
                                              grepl('111C', LRU) ~ '111C',
                                              grepl('98A', LRU) ~ '98A',
                                              TRUE ~ LRU),
                             MLRA = as.factor(MLRA)) 
lmustat <- lmustat |> group_by(mkeys, grpname,areasymbol,liid,nationalmusym,muiid,musym,muname,mukind,mutype,mustatus,dmuinvesintens,muacres,farmlndcl,dmuiid,pct_component,pct_hydric,n_component,n_majcompflag,MLRA) |> 
  summarise(n=sum(n)) |> ungroup()

mustat <- lmustat |> group_by(muiid,musym,muname,mukind,mutype,mustatus,dmuinvesintens,muacres,farmlndcl,dmuiid,pct_component,pct_hydric,n_component,n_majcompflag,MLRA) |> 
  summarise(n=sum(n)) |> ungroup()
#do the math to determine the most appropriate MLRA (MLRA affinity)
mustat <- mustat |> group_by(MLRA) |> mutate(mlratotal = sum(n), pmlra = n/mlratotal*100) |> ungroup()
mustat <- mustat |> group_by(muiid) |> mutate(mutotal = sum(n), affinity = n/mutotal*100, maxaffinity = max(affinity)) |> ungroup()
mlra_affinity <- mustat |> subset(affinity == maxaffinity)
mlra_affinity <- mlra_affinity |> group_by(muiid) |> mutate(count=length(muiid), ,
                                                            MLRAind = as.numeric(MLRA))

lmu_affinity <- mlra_affinity |> left_join(lmu[,c('lmapunitiid', 'muiid')]) |> subset(!is.na(lmapunitiid), select= c(lmapunitiid, MLRAind)) |> unique() |> as.matrix()


#Render new attributes to the mukey raster reflecting MLRA affinity
r.affinity <- classify(mkeys, lmu_affinity, others=NA)
plot(r.affinity)
writeRaster(r.affinity, 'gdbtemp/r.affinity.tif', overwrite=T)

plot(r.affinity)



# rvar <- subset(lmu,select=c('lmapunitiid','pct_hydric')) |> as.matrix()
# hydric <- classify(s.rast90, rvar, others=NA)
# plot(hydric)