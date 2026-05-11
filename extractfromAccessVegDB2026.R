library(vegnasis)
library(RODBC)
dbpath <- "C:/a/Ecological_Sites/Database/Vegetation/VegLog2026East0116.accdb"
con <- odbcConnectAccess2007(dbpath)
# sqlTables(con, tableType = "TABLE")

site <- sqlFetch(con, "Sites")
spp <- sqlFetch(con, "Observed_Species")
preveg <- site |> left_join(spp) |> group_by(Observation_ID) |> mutate(alltallshrub = sum(TallShrub)) |> ungroup() |> mutate(nstrats = ifelse(is.na(alltallshrub) | (alltallshrub == 0 & !Observation_Type %in% c('DSP', 'USFS DSP')), 4,5)) |> group_by(Observation_ID, AcTaxon) |> mutate(BA = max(BA, na.rm=T),DBH_lower = max(DBH_lower, na.rm=T),DBH_upper = max(DBH_upper, na.rm=T)) |> ungroup()
preveg <- preveg |> mutate(BA = ifelse(is.infinite(BA),NA,BA),DBH_lower = ifelse(is.infinite(DBH_lower),NA,DBH_lower),DBH_upper = ifelse(is.infinite(DBH_upper),NA,DBH_upper))

preveg1 <- preveg |> mutate(plot=ifelse(!is.na(User_Plot_ID), User_Plot_ID,Observation_ID),
                            label = Observation_Label,
                            taxon = AcTaxon,
                            cover = Field,
                            stratum.min = 0,
                            stratum.max = 0.5,
                            crown.min = Fmin,
                            crown.max = Fmax,
                            BA=BA,
                            dbh.min = NA,
                            dbh.max = NA,
                            date=(as.Date(ISOdate(Year, Mon, Day))),
                            lat=Latitude,
                            lon=Longitude) |> subset(cover > 0)
# veg <-  pre.fill.veg(preveg1)
preveg2 <- preveg |> mutate(plot=ifelse(!is.na(User_Plot_ID), User_Plot_ID,Observation_ID),
                            label = Observation_Label,
                            taxon = AcTaxon,
                            cover = Shrub,
                            stratum.min = 0.5,
                            stratum.max = ifelse(nstrats == 5, 2,5),
                            crown.min = Smin,
                            crown.max = Smax,
                            BA=BA,
                            dbh.min = NA,
                            dbh.max = NA,
                            date=(as.Date(ISOdate(Year, Mon, Day))),
                            lat=Latitude,
                            lon=Longitude) |> subset(cover > 0)
# veg <-  pre.fill.veg(preveg2)

preveg3 <- preveg |> mutate(plot=ifelse(!is.na(User_Plot_ID), User_Plot_ID,Observation_ID),
                            label = Observation_Label,
                            taxon = AcTaxon,
                            cover = TallShrub,
                            stratum.min = 2,
                            stratum.max = 5,
                            crown.min = TSmin,
                            crown.max = TSmax,
                            BA=BA,
                            dbh.min = NA,
                            dbh.max = NA,
                            date=(as.Date(ISOdate(Year, Mon, Day))),
                            lat=Latitude,
                            lon=Longitude) |> subset(cover > 0)
# veg <-  pre.fill.veg(preveg3)
preveg4 <- preveg |> mutate(plot=ifelse(!is.na(User_Plot_ID), User_Plot_ID,Observation_ID),
                            label = Observation_Label,
                            taxon = AcTaxon,
                            cover = Subcanopy,
                            stratum.min = 5,
                            stratum.max = ifelse(!is.na(Stratum8), Stratum8,15),
                            crown.min = SCmin,
                            crown.max = SCmax,
                            BA=BA,
                            dbh.min = ifelse(Tree > 0, NA,DBH_lower),
                            dbh.max = ifelse(Tree > 0, NA,DBH_upper),
                            date=(as.Date(ISOdate(Year, Mon, Day))),
                            lat=Latitude,
                            lon=Longitude) |> subset(cover > 0)
# veg <-  pre.fill.veg(preveg4)
preveg5 <- preveg |> mutate(plot=ifelse(!is.na(User_Plot_ID), User_Plot_ID,Observation_ID),
                            label = Observation_Label,
                            taxon = AcTaxon,
                            cover = Tree,
                            stratum.min = ifelse(!is.na(Stratum8) & Year > 2013, Stratum8,15),
                            stratum.max = NA,
                            crown.min = Tmin,
                            crown.max = Tmax,
                            BA=BA,
                            dbh.min = ifelse(Tree <= 0, NA,DBH_lower),
                            dbh.max = ifelse(Tree <= 0, NA,DBH_upper),
                            date=(as.Date(ISOdate(Year, Mon, Day))),
                            lat=Latitude,
                            lon=Longitude) |> subset(cover > 0)
veg <-  (rbind(preveg1,preveg2,preveg3,preveg4,preveg5)) |> subset(select=c(plot,label, Module,cover,taxon,stratum.min,stratum.max,BA,dbh.min,dbh.max,date,lat,lon)) |> unique()

#Need to average plots that use modules
anymods <- unique(veg[veg$Module %in% c('NW','NE','SW','SE'),]$plot)
vegmods <- veg |> subset(plot  %in% anymods & Module %in% c('NW','NE','SW','SE'))

vegmods <- vegmods |> group_by(plot,label, taxon,stratum.min,stratum.max,BA,dbh.min,dbh.max,date,lat,lon) |> summarise(cover = sum(cover)/4) |> ungroup() |> as.data.frame()

vegnomods <- veg |> subset(!plot  %in% anymods & Module %in% c('All'), select = -c(Module)) |> as.data.frame()
veg <- vegnomods |> rbind(vegmods[,colnames(vegnomods)])
veg <- pre.fill.veg(veg)

veg <- veg |> fill.type.df() |> fill.nativity.df() |> mutate(symbol = fill.usda.symbols(taxon)) |> fill.hts.df() 

plotids <- c('S2014MI015003', 'S2014MI081001','S2014MI117001','S2014MI123001','S2015MI015001','S2015MI037001','S2015MI045001','S2015MI045002','S2015MI107001')

veg <- subset(veg, plot %in% plotids)

veg.str <- vegnasis::get.structure(veg)
