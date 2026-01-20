library(soilDB)
library(vegnasis)
# 1	Coordinates (WGS84 lat/lon)
# 2	Geographic Elements (MLRA/State/County)
# 3	Elevation, Slope, and Landform Elements
# 4	Hydrologic Elements (drain/pond/flood)
# 5	Other Vegetation Kind Elements
# 6	Disturbance (optional)
# 7	Observers and Data Origin
# 8	Veg Plot Size
# 9	Protocol
# 10	Total Canopy Cover
# 11	Dominant Species
# 12	Cover by Sp
# 13	Strata by Sp
# 14	Complete Species Composition
# 15	Total Basal Area + by Sp
# 16	Plant Heights (max by spp + overstory base)
# 17	Basal Area Tree Diameters
# 18	Foliar Cover LPI
# 19	Ground Surface Cover
# 20	Snags
# 21	Tree Rings Counted
# 22	Tree Heights
# 23	Seedlings Counted
# 24	Biomass Documented


readcoverclass <- function(x){
  cover = case_when(
    x %in% "trace" ~ (0.1)/2,
    x %in% "0.1 to 1%" ~ (0.1+1)/2,
    x %in% "1.1 to 2%" ~ (1+2)/2,
    x %in% "2 to 5%" ~ (2+5)/2,
    x %in% "6 to 10%" ~ (5+10)/2,
    x %in% "11 to 25%" ~ (10+25)/2,
    x %in% "26 to 50%" ~ (25+50)/2,
    x %in% "51 to 75" ~ (50+75)/2,
    x %in% "76 to 95%" ~ (75+95)/2,
    x %in% "> 95%" ~ (95+100)/2,
    TRUE ~ NA)
  return(cover)}

#Load Demo Data
siteass <- vegnasis::siteass20250414
sites <- vegnasis::sites20250414
veg.raw <- vegnasis::veg.raw20250414
vegplot <- vegnasis::vegplot20250414
vegground <- vegnasis::vegground20250414

#To use your own data, remove # below
siteass <- get_site_association_from_NASIS(SS=F)
sites <- get_site_data_from_NASIS_db(SS=F)
veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS=F)
vegplot <- soilDB::get_vegplot_from_NASIS_db(SS=F)
vegground <- get_vegplot_groundsurface_from_NASIS_db(SS=F)
ba <- get_vegplot_speciesbasalarea_from_NASIS(SS=F)
si <- get_vegplot_tree_si_details_from_NASIS_db(SS=F)
trans <- get_vegplot_transect_from_NASIS_db(SS=F)
transp <- get_vegplot_transpecies_from_NASIS_db(SS=F)
siteeco <- subset(sites, select=c(usiteid, upedonid, obsdate, ecositeid, ecositenm, ecostatename, commphasename, longstddecimaldegrees,latstddecimaldegrees, horizdatnm, site_mlra, site_state, site_county,
                                  elev, slope, aspect,hillslopeprof,geomposflats,geomposhill,geompostrce, geomposmntn,
                                  drainagecl,pondfreqcl,flodfreqcl,
                                  earthcovkind1, earthcovkind2)) |> unique() 

sitepedon <- siteeco |> group_by(usiteid) |>   summarise(upedonid = first(upedonid, na_rm=T))
colnames(sitepedon) <- c('usiteid','AnyPedon')
vegplot1 <- subset(vegplot, select=c(usiteid, obsdate, assocuserpedonid, vegplotid, vegplotiid, vegplotsize, primarydatacollector,vegdataorigin, obsintensity, cancovtotalpct, cancovtotalclass, overstorycancontotalpct, overstorycancovtotalclass ,treesnagdensityhard,treesnagdensitysoft,basalareaplottotal))|> 
  mutate(overstorycover = ifelse(!is.na(overstorycancontotalpct), overstorycancontotalpct,readcoverclass(overstorycancovtotalclass)),
         allvegcover = ifelse(!is.na(cancovtotalpct), cancovtotalpct,readcoverclass(cancovtotalclass))
  )

vegground1 <- subset(vegground, select=c(vegplotid, transectlength, totalpointssampledcount, groundsurfcovtype, groundcoverptcount, groundcoverptpct)) 

ba1 <- subset(ba, select=c(vegplotid, plantsciname, speciesbasalarea, treediameterbreastheight))
trans1 <- subset(trans, select=c(vegplotid, totalpointssampledcount, groundsurfcovpointssamp, totharvestannualprod ))
transp1 <- subset(transp, select=c(vegplotid, plantsciname, speciesfoliarcovhitcount, speciesaveyielddblsamp, plantprodquadratsize))
veg <- veg.raw |> clean.veg()
checkspp <- veg |> subset(!is.na(taxon)) |> group_by(plot) |> summarize(nspp = length(taxon), covmax = max(cover, na.rm = T))
checkstrat <- veg |> subset(!is.na(stratum.min)) |> group_by(plot) |> summarize(maxstrat = max(stratum.min))
checkcrownmax <- veg |> subset(!is.na(crown.max)) |> group_by(plot) |> summarize(maxcrown = max(crown.max))
checkstructure <- veg |> fill.type.df() |> fill.hts.df() |> get.structure() |> mutate(vegcover = (1-(1-tree/100)*(1-shrub/100)*(1-herb/100))*100) 
checkba <- ba1 |> group_by(vegplotid, plantsciname) |> summarize(sppba = mean(speciesbasalarea), dbh = mean(treediameterbreastheight)) |> group_by(vegplotid) |> summarize(totalba = sum(sppba), dbh = weighted.mean(dbh, sppba)) 
checktrans <- trans1 |> group_by(vegplotid) |> summarize(npoints = sum(totalpointssampledcount),
                                                          nbiomass = mean(totharvestannualprod))
checktransp <- transp1 |> group_by(vegplotid) |> summarize(fpoints = sum(speciesfoliarcovhitcount))
checkground <- vegground1 |> group_by(vegplotid) |> summarize(gpoints = sum(groundcoverptcount))
checksi <- si |> group_by(vegplotiid) |> summarize(rings = max(growthringcount),  treeht = max(treecanopyhttop))

vegplot1 <- vegplot1 |> 
  left_join(checkspp, join_by(vegplotid==plot)) |> 
  left_join(checkstrat, join_by(vegplotid==plot)) |> 
  left_join(checkcrownmax, join_by(vegplotid==plot)) |> 
  left_join(checkstructure, join_by(vegplotid==plot)) |> 
  left_join(checkba, join_by(vegplotid==vegplotid))|> 
  left_join(checktrans, join_by(vegplotid==vegplotid)) |> 
  left_join(checktransp, join_by(vegplotid==vegplotid)) |> 
  left_join(checkground, join_by(vegplotid==vegplotid))|> 
  left_join(checksi, join_by(vegplotiid==vegplotiid))

checkplot <- siteeco |> left_join(sitepedon)|> left_join(vegplot1) 

checkplot <-  checkplot |> mutate(
  Site_ID = usiteid,
  VegPlot_ID = vegplotid, 
  Pedon_ID = case_when(is.na(assocuserpedonid) & is.na(upedonid) ~ AnyPedon,
                       is.na(assocuserpedonid) ~ upedonid,
                       TRUE ~ assocuserpedonid),
  Obs_Date = obsdate,
  Ecosite_ID = ecositeid,
  Community_Phase = commphasename,
  Obsintensity = obsintensity,
  Element1 = case_when(
    is.na(horizdatnm) & !is.na(latstddecimaldegrees) & !is.na(longstddecimaldegrees) ~ 'No Datum',
    is.na(latstddecimaldegrees) | is.na(longstddecimaldegrees) ~ 'no coordinates',
    TRUE ~ 'Pass'),
  Element2 = case_when(is.na(site_mlra)| is.na(site_state) | is.na(site_county) ~ 'Missing MLRA, County, or State',
                       TRUE ~ 'Pass'),
  Element3 = case_when(is.na(elev)| is.na(slope) | (is.na(aspect) & slope > 0) ~ 'Missing elevation, slope, or aspect',
                       (is.na(hillslopeprof) & (!is.na(geomposhill) | !is.na(geomposmntn)))  | (is.na(geomposflats) & is.na(geomposhill)& is.na(geompostrce)& is.na(geomposmntn)) ~ 'Missing one or more landform element',
                       TRUE ~ 'Pass'),
  Element4 = case_when(is.na(drainagecl) ~ 'Missing drainage class',
                       is.na(pondfreqcl) | is.na(flodfreqcl) ~ 'Missing ponding or flooding frequency',
                       TRUE ~ 'Pass'),
  Element5 = case_when(is.na(earthcovkind1) | is.na(earthcovkind2) ~ 'Missing earth cover kind',
                       TRUE ~ 'Pass'),
  Element6 = NA,

  Element7 = case_when(is.na(primarydatacollector) ~ 'Missing primary data collector',
                       is.na(vegdataorigin) ~ 'Missing data origin',
                       TRUE ~ 'Pass'),
  Element8 = case_when(is.na(vegplotsize) ~ 'Missing Veg plot size',
                       TRUE ~ 'Pass'),
  Element9 = NA,
  
  Element10 = case_when(is.na(overstorycover) ~ 'Missing overstory canopy cover',
                        is.na(allvegcover) ~ 'Missing total vegetation cover',
                        TRUE ~ 'Pass'),
  Element11 = case_when(is.na(nspp)  ~ 'Missing dominant species',
                        TRUE ~ 'Pass'),
  Element12 = case_when(is.na(covmax)  ~ 'Missing cover data',
                        TRUE ~ 'Pass'),
  Element13 = case_when(is.na(maxstrat)  ~ 'Missing strata',
                        TRUE ~ 'Pass'),
  Element14 = case_when((vegcover < 0.5*allvegcover) | (tree < 0.5*overstorycover) ~ 'Likely incomplete species list',
                        TRUE ~ 'Pass'),
  Element15 = case_when(is.na(basalareaplottotal) ~ 'Missing total basal area',
                        is.na(totalba) ~ 'Missing species basal area',
                        abs(totalba - basalareaplottotal) > 0 ~ 'Inconsistent basal area',
                        TRUE ~ 'Pass'),
  Element16 = case_when(is.na(maxcrown) ~ 'Missing maximum vegetation height',
                        (maxcrown > 60 & !site_state %in% c('CA','OR','WA'))| maxcrown > 115  ~ 'Maximum vegetation height likely too high',
                        maxcrown <= maxstrat ~ 'Maximum vegetation height too low',
                        abs(totalba - basalareaplottotal) > 0 ~ 'Inconsistent basal area',
                        TRUE ~ 'Pass'),
  Element17 = case_when(is.na(dbh) ~ 'Missing tree diameters',
                        TRUE ~ 'Pass'),
  Element18 = case_when(is.na(npoints) ~ 'Missing total number of points',
                        is.na(fpoints) ~ 'Missing foliar cover counts',
                        TRUE ~ 'Pass'),
  Element19 = case_when(is.na(npoints) ~ 'Missing total number of points',
                        is.na(gpoints) ~ 'Missing ground surface cover counts',
                        TRUE ~ 'Pass'),
  Element20 = case_when(is.na(treesnagdensityhard) & is.na(treesnagdensitysoft) ~ 'Missing snag count',
                        TRUE ~ 'Pass'),
  Element21 = case_when(is.na(rings) ~ 'Missing ring count',
                        rings > 4500 ~ 'Likely too old',
                        TRUE ~ 'Pass'),
  Element22 = case_when(is.na(treeht) ~ 'Missing tree heights',
                        (treeht > 60 & !site_state %in% c('CA','OR','WA'))| treeht > 115  ~ 'Maximum tree height likely too high',
                        TRUE ~ 'Pass'),
  Element23 = NA,
  ,
  Element24 = case_when(is.na(nbiomass)  ~ 'Missing biomass',
                        TRUE ~ 'Pass')
) |> subset(!is.na(VegPlot_ID), select=c("Site_ID", "VegPlot_ID","Pedon_ID","Obs_Date",
                     "Ecosite_ID","Community_Phase","Obsintensity","Element1",
                     "Element2","Element3","Element4","Element5",
                     "Element6","Element7","Element8","Element9",
                     "Element10","Element11","Element12","Element13",
                     "Element14","Element15","Element16","Element17",
                     "Element18","Element19","Element20","Element21",
                     "Element22","Element23","Element24"))

rm(siteass,siteeco,checkstrat,checkspp,checkcrownmax,checkstructure,checkba,checktrans,checktransp,checkground,checksi)
write.csv(checkplot, 'checkplot.csv', row.names = F, na="")

checkplot2 <- checkplot |> mutate(Element1 =  ifelse(Element1 %in% "Pass", 1,ifelse(is.na(Element1),NA,0)),
                                  Element2 =  ifelse(Element2 %in% "Pass", 1,ifelse(is.na(Element2),NA,0)),
                                  Element3 =  ifelse(Element3 %in% "Pass", 1,ifelse(is.na(Element3),NA,0)),
                                  Element4 =  ifelse(Element4 %in% "Pass", 1,ifelse(is.na(Element4),NA,0)),
                                  Element5 =  ifelse(Element5 %in% "Pass", 1,ifelse(is.na(Element5),NA,0)),
                                  Element6 =  ifelse(Element6 %in% "Pass", 1,ifelse(is.na(Element6),NA,0)),
                                  Element7 =  ifelse(Element7 %in% "Pass", 1,ifelse(is.na(Element7),NA,0)),
                                  Element8 =  ifelse(Element8 %in% "Pass", 1,ifelse(is.na(Element8),NA,0)),
                                  Element9 =  ifelse(Element9 %in% "Pass", 1,ifelse(is.na(Element9),NA,0)),
                                  Element10 =  ifelse(Element10 %in% "Pass", 1,ifelse(is.na(Element10),NA,0)),
                                  Element11 =  ifelse(Element11 %in% "Pass", 1,ifelse(is.na(Element11),NA,0)),
                                  Element12 =  ifelse(Element12 %in% "Pass", 1,ifelse(is.na(Element12),NA,0)),
                                  Element13 =  ifelse(Element13 %in% "Pass", 1,ifelse(is.na(Element13),NA,0)),
                                  Element14 =  ifelse(Element14 %in% "Pass", 1,ifelse(is.na(Element14),NA,0)),
                                  Element15 =  ifelse(Element15 %in% "Pass", 1,ifelse(is.na(Element15),NA,0)),
                                  Element16 =  ifelse(Element16 %in% "Pass", 1,ifelse(is.na(Element16),NA,0)),
                                  Element17 =  ifelse(Element17 %in% "Pass", 1,ifelse(is.na(Element17),NA,0)),
                                  Element18 =  ifelse(Element18 %in% "Pass", 1,ifelse(is.na(Element18),NA,0)),
                                  Element19 =  ifelse(Element19 %in% "Pass", 1,ifelse(is.na(Element19),NA,0)),
                                  Element20 =  ifelse(Element20 %in% "Pass", 1,ifelse(is.na(Element20),NA,0)),
                                  Element21 =  ifelse(Element21 %in% "Pass", 1,ifelse(is.na(Element21),NA,0)),
                                  Element22 =  ifelse(Element22 %in% "Pass", 1,ifelse(is.na(Element22),NA,0)),
                                  Element23 =  ifelse(Element23 %in% "Pass", 1,ifelse(is.na(Element23),NA,0)),
                                  Element24 =  ifelse(Element24 %in% "Pass", 1,ifelse(is.na(Element24),NA,0)))
colnames(checkplot2) <- c("Site_ID", "VegPlot_ID","Pedon_ID","Obs_Date",
                          "Ecosite_ID","Community_Phase","Obsintensity",
                          "Coordinates (WGS84 lat/lon)",
                                     "Geographic Elements (MLRA/State/County)",
                                     "Elevation, Slope, and Landform Elements",
                                     "Hydrologic Elements (drain/pond/flood)",
                                     "Other Vegetation Kind Elements",
                                     "Disturbance (optional)",
                                     "Observers and Data Origin",
                                     "Veg Plot Size",
                                     "Protocol",
                                     "Total Canopy Cover",
                                     "Dominant Species",
                                     "Cover by Sp",
                                     "Strata by Sp",
                                     "Complete Species Composition",
                                     "Total Basal Area + by Sp",
                                     "Plant Heights (max by spp + overstory base)",
                                     "Basal Area Tree Diameters",
                                     "Foliar Cover LPI",
                                     "Ground Surface Cover",
                                     "Snags",
                                     "Tree Rings Counted",
                                     "Tree Heights",
                                     "Seedlings Counted",
                                     "Biomass Documented")

write.csv(checkplot2, 'checkplot2.csv', row.names = F, na="")
