library(soilDB)
library(vegnasis)
library(sf)
library(xlsx)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts <- st_read('SoilSiteEcotypeVeg.shp')
pts <- st_transform(pts, crs='EPSG:4326')
pts <- pts |>  mutate(lat=st_coordinates(pts)[,2],lon=st_coordinates(pts)[,1]) |> st_drop_geometry()
plnts <- read.csv('SoilSiteEcotypeVeg_PlantsData.csv')
sls <- read.csv('SoilSiteEcotypeVeg_SoilsData.csv')
chek <- read.csv('redoout.csv')
siteids <- read.csv('site_ids.csv')
siteids <- unique(siteids[,c('NASIS.SiteID','Project_SiteEcotypeVegSoil.recordid','Lat','Long')]) |> subset(!NASIS.SiteID %in% 'SKIP' & !is.na(NASIS.SiteID))
#siteids <- siteids |> group_by(NASIS.SiteID) |> mutate(nrecs = length(NASIS.SiteID))
pts <- pts |> left_join(siteids, join_by(recordid==Project_SiteEcotypeVegSoil.recordid)) |> subset(!NASIS.SiteID %in% 'SKIP' & !is.na(NASIS.SiteID) & !is.na(datecollec))
ptscomment <- pts |> subset(select = c(NASIS.SiteID, commentgen))


strats1 <- plnts |> subset(!is.na(overstory) & nchar(overstory) > 0) |> mutate(stratum.min = 5, stratum.max = NA, cover = overstory)
strats2 <- plnts |> subset(!is.na(understory) & nchar(understory) > 0) |> mutate(stratum.min = 0.5, stratum.max = 5, cover = understory)
strats3 <- plnts |> subset(!is.na(groundcover) & nchar(groundcover) > 0) |> mutate(stratum.min = 0, stratum.max = 0.5, cover = groundcover)
strats <- rbind(strats1,strats2,strats3)
strats <- strats |> mutate(cover = case_when(cover %in% 'A' ~ 50, #>20%
                                             cover %in% 'C' ~ 10, #2-20%
                                             cover %in% 'R' ~ 1)) #<2%
rm(strats1,strats2,strats3)
pts <- unique(pts[,c('NASIS.SiteID','recordid','datecollec','lat','lon','observer')])
#pts <- pts |> group_by(NASIS.SiteID) |> mutate(nrecs = length(NASIS.SiteID))

strats <- strats |> left_join(pts, join_by(record_id==recordid))
strats <- strats |> left_join(chek[,c('misname','final')], join_by(plants_code==misname))
strats <- strats |> mutate(symbol = ifelse(is.na(final), plants_code, final), taxon = vegnasis::fill.taxon.from.symbols(symbol), chek = check.phytogeography(taxon, 'MI'), plot=record_id, date=datecollec, label = plants_code)

veg <- strats |> vegnasis::pre.fill.veg()
veg <- veg |> mutate(taxon = harmonize.taxa(veg$taxon, fix = TRUE, sensu = "usda")) |> fill.type.df() |> fill.nativity.df() |> fill.hts.df()
veg <- veg |> mutate(symbol = ifelse(!type %in% c('tree','shrub/vine') & stratum.min > 0.5, '2TREE', symbol),
                     taxon = ifelse(!type %in% c('tree','shrub/vine') & stratum.min > 0.5, 'Tree', taxon),
                     type = ifelse(!type %in% c('tree','shrub/vine') & stratum.min > 0.5, 'tree', type))
veg <- veg |> left_join(unique(pts[!is.na(pts$observer),c('recordid','observer')]), join_by(plot==recordid)) |> left_join(unique(pts[!is.na(pts$observer),c('recordid','NASIS.SiteID')]), join_by(plot==recordid)) 
veg <- veg |> mutate(plot = NASIS.SiteID) |> subset(!is.na(plot), select=-c(NASIS.SiteID))
veg <- veg |> mutate(siteid=plot, vegplotid=plot,vegplotname="", plotsize = NA, obsintensity = "low-intensity",
                     #obskind = "actual site observation date",dataorigin="spreadsheet form",
                     treecover = NA, vegcover=NA, totalBA=NA, ES.ID=NA, ES.state=NA, ES.phase=NA)
veg <- subset(veg, select=c("siteid","lat","lon","date","vegplotid","vegplotname","observer","plotsize", "obsintensity", "treecover","vegcover", "totalBA", 
                            "symbol", "taxon","type","nativity","cover","stratum.min", "stratum.max", "crown.min","crown.max","dbh.min","dbh.max","BA")) |> arrange(vegplotid, symbol, stratum.min) |> mutate(stratum.min = vegnasis::ht.USC(stratum.min), stratum.max = vegnasis::ht.USC(stratum.max))

theseplots <- veg$vegplotid |> unique()
for(i in 1:length(theseplots)){
thisplot <- theseplots[i]

theserows <- subset(veg, vegplotid %in% thisplot)

# Create an empty workbook
wb <- createWorkbook()

# Create sheets within the workbook
exportNASIS <- createSheet(wb, sheetName = "exportNASIS")

# Add OSD data
addDataFrame(
  theserows,
  exportNASIS,
  startRow = 3,
  startColumn = 1,
  row.names = F
)

# Create a data.frame with the name and version of the workbook required for NASIS
wbnamecell <- data.frame("Simple Releves", "1.0")

# Add the required NASIS data.frame to the workbook
addDataFrame(
  wbnamecell,
  exportNASIS,
  startRow = 1,
  startColumn = 1,
  col.names = F,
  row.names = F
)

# Save the workbook
saveWorkbook(wb, paste0("exports/vegplot-", thisplot, ".xlsx"))
}














write.csv(veg, 'exportNASIS.csv', na='', row.names = F)
write.csv(ptscomment, 'ptscomment.csv', na='', row.names = F)


ptscomment <- ptscomment |> mutate(maxht = stringr::str_extract(commentgen))
























# #QC plant symbols
# strats <- strats |> mutate(taxon = vegnasis::fill.taxon.from.symbols(plants_code), chek = check.phytogeography(taxon, 'MI'))
# missing <-  subset(strats, !chek %in% 'pass')
# redo <- data.frame(misname = unique(missing$plants_code))
# redo$alt = NA
# redo$stdist = NA
# usda <-  vegnasis::usdaplants |> mutate(chek = check.phytogeography(taxon, 'MI')) |> subset(chek %in% 'pass')
# for(i in 1:nrow(redo)){#i=2
# thissp <- redo$misname[i]
# usda1 <- usda |> mutate(std = stringdist::stringdist(thissp, sym, method = 'jw'), mindist = min(std)) |> subset(mindist == std)
# redo$alt[i] <- paste(usda1$sym,  collapse = '; ')
# redo$stdist = min(usda1$std)}
# redo$alt <- substr(redo$alt, 1,200)
# write.csv(redo, 'redo.csv',row.names = F)
other <- vegnasis::usdaplants

usda <-  vegnasis::usdaplants |> mutate(chek = check.phytogeography(taxon, 'MI')) |> subset(chek %in% 'pass')

thissp <- 'DIPO'
usda1 <- usda |> mutate(std = stringdist::stringdist(thissp, sym, method = 'jw'), mindist = min(std))
theseplots <- strats |> subset(record_id %in% strats[strats$plants_code %in% thissp,]$record_id)
