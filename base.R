library(soilDB)
library(stringr)
library(dplyr)

#install and load package
remotes::install_github("phytoclast/vegnasis", dependencies = FALSE)
library(vegnasis)
#fresh NASIS data
veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS = FALSE)
veg.raw2 <- soilDB::get_vegplot_transpecies_from_NASIS_db(SS = FALSE)
site.raw <- soilDB::get_site_data_from_NASIS_db(SS = FALSE)
saveRDS(veg.raw, 'data/nasisvegraw2.RDS')
saveRDS(veg.raw2, 'data/nasisvegtransraw2.RDS')
saveRDS(site.raw, 'data/nasissiteraw2.RDS')


veg.raw3 <- soilDB::get_vegplot_transect_from_NASIS_db (SS = FALSE)
veg <- clean.veg.transect(veg.raw2) |> fill.hts.df()
#sample data
veg.raw <- vegnasis::nasis.veg
veg <- clean.veg(veg.raw)
#fill in missing data
veg <- clean.veg(veg.raw) |> fill.type.df()  |> fill.hts.df() |> fill.nativity.df()
#get wetland status
hydric <- get.wetness(veg, region = 'NCNE')
#get structure and dominance association
structure <- get.structure(veg)
association <- get.assoc(veg)
#quantified structural breakdown based on stratum identity and using USNVC defaults
strat.summary <- vegnasis::summary.strata(veg)
#quantified structural breakdown based on live crown and using EDIT's breaks
breaks.ft <- c(0.5, 1, 2, 4.5, 13, 40, 80, 120)
breaks <- vegnasis::ht.metric(breaks.ft)
EDIT.strat.summary <- vegnasis::summary.crown.thickness(veg, breaks)

#prepare for analysis
m <- make.plot.matrix(veg, tr = 'log')
d = vegan::vegdist(m, method='bray')
t <- cluster::agnes(d, method = 'ward')|> as.hclust()
ape::plot.phylo(ape::as.phylo(t))
t <- optpart::flexbeta(d, beta = -0.2)|> as.hclust()
ape::plot.phylo(ape::as.phylo(t))


veg$new <- get.habit.code(veg$plantsciname)
veg$new2 <- get.habit(veg$new)
m <- make.plot.matrix(veg)
# remotes::install_github("natearoe/ecositer", dependencies = FALSE)
# library(ecositer)
# test <- veg_summary(veg_df = ecositer::vegetation_dataframe)
#
# test.R022AB006CA <- test$R022AB006CA$Raw_data

veg.raw <- readRDS('data_raw/veg.raw.select.RDS')
usethis::use_data(veg.raw, overwrite = T)
# taxon.habits <- read.csv('data_raw/taxon.habits.csv')
# usethis::use_data(taxon.habits, overwrite = T)
# gho <- read.csv('data_raw/gho.csv')
# usethis::use_data(gho, overwrite = T)
# genus.habits <- read.csv('data_raw/genus.habits.csv')
# usethis::use_data(genus.habits, overwrite = T)
# syns <- read.csv('data_raw/m.ac.csv')
# usethis::use_data(syns, overwrite = T)
# nasis.veg <- readRDS('data_raw/veg.spp.RDS')
# usethis::use_data(nasis.veg, overwrite = T)
# hydric <- read.csv('data_raw/hydric.csv')
# usethis::use_data(hydric, overwrite = T)
# usdaplants <- read.csv('data_raw/plantssym.csv')
# PLANTS <- read.csv('data_raw/PLANTSdownloadData.txt')
# PLANTS.illegit <- PLANTS |> subset(grepl('auct.',Genera.Binomial.Author) |
#                                      grepl('illeg.',Genera.Binomial.Author) |
#                                      grepl(' non',Genera.Binomial.Author) |
#                                      grepl('auct',Trinomial.Author) |
#                                      grepl('illeg.',Trinomial.Author) |
#                                      grepl(' non',Trinomial.Author),
#                                    select=c(Accepted.Symbol, Symbol, Scientific.Name, Genera.Binomial.Author, Trinomial.Author))
# usdaplants <- usdaplants |> subset(!sym %in% PLANTS.illegit$Symbol)
# usethis::use_data(usdaplants, overwrite = T)
# natdat <- read.csv('data_raw/nativity.csv')
# usethis::use_data(natdat, overwrite = T)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

veg.plot <- readRDS('data/veg.plot.RDS')
veg.data <-  readRDS('data/veg.data.RDS')
veg.loc <-  readRDS('data/veg.loc.RDS')
veg.spp <-  readRDS('data/veg.spp.RDS')
veg.data.veg <-  readRDS('data/veg.data.veg.RDS')
veg.data.trans <-  readRDS('data/veg.data.trans.RDS')

plant.hts <- read.delim('data/plants/Plant_heights.txt')
veg <- clean.veg(veg.spp)
veg <- fill.hts.df(veg)
breaks <- c(0.1, 0.5, 2, 5, 10, 20, 30)
strat.summary <- vegnasis::summary.crown.thickness(veg, breaks)

#
# veg.plot <- soilDB::get_vegplot_from_NASIS_db()
# veg.data <- soilDB::get_veg_data_from_NASIS_db()
# veg.loc <- soilDB::get_vegplot_location_from_NASIS_db()
# veg.spp <- soilDB::get_vegplot_species_from_NASIS_db()
# veg.data.veg <- veg.data$veg
# veg.data.trans <- veg.data$vegtransect
#
# saveRDS(veg.plot, 'data/veg.plot.RDS')
# saveRDS(veg.data, 'data/veg.data.RDS')
# saveRDS(veg.loc, 'data/veg.loc.RDS')
# saveRDS(veg.spp, 'data/veg.spp.RDS')
# saveRDS(veg.data.veg, 'data/veg.data.veg.RDS')
# saveRDS(veg.data.trans, 'data/veg.data.trans.RDS')

ht.round <- function(ht){
         ifelse(ht >= 8,round(ht,0),
                ifelse(ht >= 3, floor(ht*2+0.499)/2,round(ht,1)
                ))}


ht.metric <- function(ft){
  round(as.numeric(ft)*0.3048,1)
}
ht.medieval <- function(m){
  round(as.numeric(m)/0.3048,1)
}
diam.metric <- function(inch){
  round(as.numeric(inch)*2.54,0)
}
diam.medieval <- function(cm){
  round(as.numeric(cm)/2.54,1)
}

syns <- read.csv('data/plants/m.ac.csv')

cover.agg <- function(x){round(100*(1-10^(sum(log10(1-(x/100.001))))),1)}

clean.veg <- function(x){
  plant.hts <- read.delim('data/plants/Plant_heights.txt')

  # x <- subset(x, select = c(vegplotid, plantsym, plantsciname, planttypegroup, plantnativity, plantheightcllowerlimit, plantheightclupperlimit, livecanopyhtbottom, livecanopyhttop,
  #                                 overstorydbhmin, overstorydbhmax, speciestraceamtflag,
  #                                 speciescancovpct, speciescancovclass, speciescomppct, speciescompbywtpct,
  #                                 akstratumcoverclass, akstratumcoverclasspct))

  x <- x %>% left_join(plant.hts, by=c('plantsciname'='Scientific.Name'))
  x <- x %>% mutate(Ht_m = case_when(
    !is.na(Ht_m) ~ Ht_m,
    planttypegroup %in%  "tree" ~ 25,
    planttypegroup %in%  "shrub/vine" ~ 3,
    planttypegroup %in%  "forb" ~ 0.6,
    planttypegroup %in%  "grass/grasslike" ~ 0.6,
    planttypegroup %in%  "moss" ~ 0,
    TRUE ~ 0))
  x <- x %>% mutate(
    cover = case_when(
      !is.na(akstratumcoverclasspct) ~ as.numeric(akstratumcoverclasspct),
      !is.na(speciescancovpct) ~ as.numeric(speciescancovpct) + ifelse(speciestraceamtflag,0.2,0),
      speciescancovclass %in% "trace" ~ (0.1)/2,
      speciescancovclass %in% "0.1 to 1%" ~ (0.1+1)/2,
      speciescancovclass %in% "1.1 to 2%" ~ (1+2)/2,
      speciescancovclass %in% "2 to 5%" ~ (2+5)/2,
      speciescancovclass %in% "6 to 10%" ~ (5+10)/2,
      speciescancovclass %in% "11 to 25%" ~ (10+25)/2,
      speciescancovclass %in% "26 to 50%" ~ (25+50)/2,
      speciescancovclass %in% "51 to 75" ~ (50+75)/2,
      speciescancovclass %in% "76 to 95%" ~ (75+95)/2,
      speciescancovclass %in% "> 95%" ~ (95+100)/2,
      !is.na(speciescomppct) ~ as.numeric(speciescomppct),
      TRUE ~ 0),

    ht.min = 0,

    ht.max = case_when(
      !is.na(livecanopyhttop) ~ ht.metric(livecanopyhttop),
      !is.na(plantheightclupperlimit) ~ pmin(ht.metric(plantheightclupperlimit),
                                             pmax(Ht_m, ht.metric(plantheightcllowerlimit) +
                                                    pmax(Ht_m, ht.metric(plantheightcllowerlimit)/5))),
      akstratumcoverclass %in% "tree regeneration generally less than 4.5 m (15 ft) tall" ~ 4.5,
      akstratumcoverclass %in% "stunted tree generally less than 4.5 m (15 ft) tall" ~ 4.5,
      akstratumcoverclass %in% "medium tree generally between 4.5 and 12 m (15 and 40 ft) tall" ~ 12,
      akstratumcoverclass %in% "tall tree generally greater than 12 m (40 ft) tall" ~
        pmax(Ht_m*0.9, 12+3),
      akstratumcoverclass %in% "dwarf shrub layer less than about 20 cm (8 in) tall" ~ 0.2,
      akstratumcoverclass %in% "low shrub between about 20 and 100 cm (8 and 36 in) tall" ~ 1,
      akstratumcoverclass %in% "medium shrub between about 1 and 3 m (3 and 10 ft) tall" ~ 3,
      akstratumcoverclass %in% "tall shrub greater than about 3 m (10 ft) tall" ~
        pmin(5,pmax(Ht_m*0.9, 3+1)),
      akstratumcoverclass %in% "low and dwarf graminoid less than about 10 cm (4 in) tall" ~ 0.1,
      akstratumcoverclass %in% "medium graminoid between about 10 and 60 cm (4 and 24 in) tall" ~ 0.6,
      akstratumcoverclass %in% "tall graminoid generally greater than 60 cm (24 in) tall" ~
        pmin(2,pmax(Ht_m*0.9, 0.6+0.6)),
      akstratumcoverclass %in% "low and dwarf forb generally less than 10 cm (4 in) tall" ~ 0.1,
      akstratumcoverclass %in% "medium forb between about 10 and 60 cm (4 and 24 in) tall" ~ 0.6,
      akstratumcoverclass %in% "tall forb generally greater than 60 cm (24 in) tall" ~
        pmin(2,pmax(Ht_m*0.9, 0.6+0.6)),
      akstratumcoverclass %in% "mosses" ~ 0,
      TRUE ~ NA_real_),

    ht.min = case_when(
      !is.na(livecanopyhtbottom) ~ ht.metric(livecanopyhtbottom),
      TRUE ~ ht.max/2),

    ht.max = ht.round(ht.max),
    ht.min = ht.round(ht.min),
    diam.min = diam.metric(overstorydbhmin),
    diam.max = diam.metric(overstorydbhmax))
  colnames(x)
  x <- x %>% subset(select= c("vegplotid","plantsym","plantsciname","planttypegroup",
              "plantnativity","cover","ht.min","ht.max","diam.min","diam.max"))
  return(x)
}



x <- clean.veg(veg.spp)

y<- x %>% group_by(vegplotid, planttypegroup) %>% summarise(Cover = cover.agg(cover))

#stratum summary ----
# breaks <- c(0.1, 0.5, 2, 5, 10, 20, 30)
# summary.strata <-  function(x, breaks){
#   x <- x %>% mutate(stratum=0)
#   for(i in 1:length(breaks)){
#     x <- x %>% mutate(stratum= ifelse(ht.max >= breaks[i],breaks[i], stratum))
#   }
#   y <- x %>% group_by(vegplotid, planttypegroup, stratum) %>% summarise(Cover = cover.agg(cover))
#   return(y)
# }
# y <- summary.strata(x, breaks)


#stratum summary ----
breaks <- c(0.1, 0.5, 2, 5, 10, 20, 30)
summary.strata <-  function(x, breaks){
  nbks <- length(breaks)+1
  brks <- c(0,breaks,1000)
  for(i in 1:(nbks)){#i = 8
    y0 <- x %>% subset(ht.max < brks[i+1] & ht.max >= brks[i])

    if(nrow(y0)>0){
      y0 <- y0 %>% mutate(stratum=i, stratum.label = paste0(brks[i], "-", ifelse(i==nbks, "+",brks[i+1])))
      y1 <- y0 %>% group_by(vegplotid, planttypegroup, stratum, stratum.label) %>% summarise(Cover = cover.agg(cover))
    }
    if(i==1){y <- y1}else{y <- rbind(y, y1)}
  }
  return(y)
}

y1 <- summary.strata(x, breaks)

#live crown thickness ----
breaks <- c(0.1, 0.5, 2, 5, 10, 20, 30)
summary.crown.thickness <-  function(x, breaks){
  nbks <- length(breaks)+1
  brks <- c(0,breaks,1000)
  for(i in 1:(nbks)){#i = 5
    y0 <- x %>% subset(ht.min < brks[i+1] & ht.max >= brks[i])

    if(nrow(y0)>0){
      y0 <- y0 %>% mutate(stratum=i, stratum.label = paste0(brks[i], "-", ifelse(i==nbks, "+",brks[i+1])))
      y1 <- y0 %>% group_by(vegplotid, planttypegroup, stratum, stratum.label) %>% summarise(Cover = cover.agg(cover))
    }
    if(i==1){y <- y1}else{y <- rbind(y, y1)}
  }
  return(y)
}

y2 <- summary.crown.thickness(x, breaks)


x <- veg.spp


clean.veg <- function(x){
  x <- x %>% mutate(
    cover = case_when(
      !is.na(akstratumcoverclasspct) ~ as.numeric(akstratumcoverclasspct),
      !is.na(speciescancovpct) ~ as.numeric(speciescancovpct) + ifelse(speciestraceamtflag,0.2,0),
      speciescancovclass %in% "trace" ~ (0.1)/2,
      speciescancovclass %in% "0.1 to 1%" ~ (0.1+1)/2,
      speciescancovclass %in% "1.1 to 2%" ~ (1+2)/2,
      speciescancovclass %in% "2 to 5%" ~ (2+5)/2,
      speciescancovclass %in% "6 to 10%" ~ (5+10)/2,
      speciescancovclass %in% "11 to 25%" ~ (10+25)/2,
      speciescancovclass %in% "26 to 50%" ~ (25+50)/2,
      speciescancovclass %in% "51 to 75" ~ (50+75)/2,
      speciescancovclass %in% "76 to 95%" ~ (75+95)/2,
      speciescancovclass %in% "> 95%" ~ (95+100)/2,
      !is.na(speciescomppct) ~ as.numeric(speciescomppct),
      TRUE ~ 0),


    strat.ht.max = case_when(
      !is.na(plantheightclupperlimit) ~ ht.metric(plantheightclupperlimit),
      akstratumcoverclass %in% "tree regeneration generally less than 4.5 m (15 ft) tall" ~ 4.5,
      akstratumcoverclass %in% "stunted tree generally less than 4.5 m (15 ft) tall" ~ 4.5,
      akstratumcoverclass %in% "medium tree generally between 4.5 and 12 m (15 and 40 ft) tall" ~ 12,
      akstratumcoverclass %in% "tall tree generally greater than 12 m (40 ft) tall" ~ NA_real_,
      akstratumcoverclass %in% "dwarf shrub layer less than about 20 cm (8 in) tall" ~ 0.2,
      akstratumcoverclass %in% "low shrub between about 20 and 100 cm (8 and 36 in) tall" ~ 1,
      akstratumcoverclass %in% "medium shrub between about 1 and 3 m (3 and 10 ft) tall" ~ 3,
      akstratumcoverclass %in% "tall shrub greater than about 3 m (10 ft) tall" ~ NA_real_,
      akstratumcoverclass %in% "low and dwarf graminoid less than about 10 cm (4 in) tall" ~ 0.1,
      akstratumcoverclass %in% "medium graminoid between about 10 and 60 cm (4 and 24 in) tall" ~ 0.6,
      akstratumcoverclass %in% "tall graminoid generally greater than 60 cm (24 in) tall" ~ NA_real_,
      akstratumcoverclass %in% "low and dwarf forb generally less than 10 cm (4 in) tall" ~ 0.1,
      akstratumcoverclass %in% "medium forb between about 10 and 60 cm (4 and 24 in) tall" ~ 0.6,
      akstratumcoverclass %in% "tall forb generally greater than 60 cm (24 in) tall" ~ NA_real_,
      akstratumcoverclass %in% "mosses" ~ 0,
      TRUE ~ NA_real_),

    strat.ht.min = case_when(
      !is.na(plantheightcllowerlimit) ~ ht.metric(plantheightcllowerlimit),
      akstratumcoverclass %in% "tree regeneration generally less than 4.5 m (15 ft) tall" ~ 0,
      akstratumcoverclass %in% "stunted tree generally less than 4.5 m (15 ft) tall" ~ 0,
      akstratumcoverclass %in% "medium tree generally between 4.5 and 12 m (15 and 40 ft) tall" ~ 4.5,
      akstratumcoverclass %in% "tall tree generally greater than 12 m (40 ft) tall" ~ 12,
      akstratumcoverclass %in% "dwarf shrub layer less than about 20 cm (8 in) tall" ~ 0,
      akstratumcoverclass %in% "low shrub between about 20 and 100 cm (8 and 36 in) tall" ~ 0.2,
      akstratumcoverclass %in% "medium shrub between about 1 and 3 m (3 and 10 ft) tall" ~ 1,
      akstratumcoverclass %in% "tall shrub greater than about 3 m (10 ft) tall" ~ 3,
      akstratumcoverclass %in% "low and dwarf graminoid less than about 10 cm (4 in) tall" ~ 0,
      akstratumcoverclass %in% "medium graminoid between about 10 and 60 cm (4 and 24 in) tall" ~ 0.1,
      akstratumcoverclass %in% "tall graminoid generally greater than 60 cm (24 in) tall" ~ 0.6,
      akstratumcoverclass %in% "low and dwarf forb generally less than 10 cm (4 in) tall" ~ 0,
      akstratumcoverclass %in% "medium forb between about 10 and 60 cm (4 and 24 in) tall" ~ 0.1,
      akstratumcoverclass %in% "tall forb generally greater than 60 cm (24 in) tall" ~ 0.6,
      akstratumcoverclass %in% "mosses" ~ 0,
      TRUE ~ NA_real_),

    live.ht.min = ht.metric(livecanopyhtbottom),
    live.ht.max = ht.metric(livecanopyhttop),

    diam.min = diam.metric(overstorydbhmin),
    diam.max = diam.metric(overstorydbhmax))

  x <- x %>% subset(select= c("vegplotid","plantsym","plantsciname","planttypegroup",
                              "plantnativity","cover","strat.ht.min","strat.ht.max","live.ht.min","live.ht.max","diam.min","diam.max"))
  return(x)
}





veg <- clean.veg(veg.spp)
veg2 <- fill.hts.df(veg)
veg3 <- vegnasis::summary.crown.thickness(veg2, c(0.5,2,5,10))
x <- veg
plot <- x$vegplotid
taxon <- x$taxon
type <- x$type
stratum.min <- x$stratum.min
stratum.max <- x$stratum.max
crown.min <- x$crown.min
crown.max <- x$crown.max

harmonize.heights <- function(plot = NA_character_ ,
                              taxon = NA_character_ ,
                              type = NA_character_ ,
                              stratum.min = NA_real_,
                              stratum.max = NA_real_,
                              crown.min = NA_real_,
                              crown.max = NA_real_){
taxon.max = get.ht.max(taxon)

df <- data.frame(plot, taxon, type, stratum.min, stratum.max, crown.min, crown.max, taxon.max) |> as.data.frame()

df <- df |> group_by(plot) |> mutate(stand.max = pmax(max(stratum.max, na.rm = TRUE, warnings =FALSE), max(crown.max, na.rm = TRUE, warnings =FALSE), na.rm = TRUE, warnings =FALSE),
                                     base.max = pmax(max(stratum.min, na.rm = TRUE, warnings =FALSE), max(crown.min, na.rm = TRUE, warnings =FALSE), na.rm = TRUE, warnings =FALSE ),
                                     stand.max =  ifelse(is.na(stand.max) | base.max >= stand.max, NA_real_, stand.max),
                                     base.max = NULL)
df <- df |> mutate(ht.max = case_when(
  !is.na(crown.max) ~ crown.max,
  !is.na(stratum.max) ~ pmin(stratum.max, taxon.max),
  TRUE ~ pmin(taxon.max, stand.max)),
  ht.max = case_when(
    !is.na(ht.max) ~ ht.max,
    type %in%  "tree" ~ 24,
    type %in%  "shrub/vine" ~ 3,
    type %in%  "forb" ~ 0.6,
    type %in%  "grass/grasslike" ~ 0.6,
    type %in%  "moss" ~ 0,
    TRUE ~ 0))

df <- df |> mutate(
  ht.max = ifelse(!is.na(stratum.min) & stratum.min >=  ht.max, stratum.min + (stratum.max-stratum.min)/10, ht.max),
                   ht.min =  case_when(
                     !is.na(crown.min) ~ crown.min,
                     TRUE ~ ht.max/2))
return(data.frame(ht.min = ht.round(df$ht.min),ht.max = ht.round(df$ht.max)))
}


hts <- harmonize.heights(plot = plot ,
                  taxon = taxon ,
                  type = type ,
                  stratum.min = stratum.min,
                  stratum.max = stratum.max,
                  crown.min = crown.min,
                  crown.max = crown.max)

system.time(
x <- veg |> mutate(stratum = case_when(
  (ht.max > 15 & (type %in% c('tree') | is.na(type)))| (type %in% 'tree' & is.na(ht.max)) ~ 'tree2',
  (ht.max > 5 & (type %in% c('shrub/vine', 'tree') | is.na(type)))| (type %in% 'tree' & is.na(ht.max)) ~ 'tree1',
  (ht.max > 0.5 & (type %in% c('shrub/vine', 'tree') | is.na(type)))| (type %in% 'shrub/vine' & is.na(ht.max)) ~ 'shrub',
  (ht.max > 0 & (type %in% c('shrub/vine') | is.na(type)))| (type %in% c('forb','grass/grasslike')) ~ 'field',
  !type %in% c('tree', 'shrub/vine', 'forb','grass/grasslike') ~ 'ground',
  TRUE ~ 'none'))
)
system.time(
  x <- veg |> mutate(stratum = case_when(
    is.na(ht.max) | is.na(type) ~ 'excluded',
    ht.max > 15 & type %in% c('shrub/vine', 'tree')  ~ '4',
    ht.max > 5  & type %in% c('shrub/vine', 'tree')  ~ '3',
    ht.max > 0.5 & type %in% c('shrub/vine', 'tree') ~ '2',
    ht.max > 0 & type %in% c('shrub/vine', 'forb','grass/grasslike') ~ '1',
    !type %in% c('tree', 'shrub/vine', 'forb','grass/grasslike') ~ '0',
    TRUE ~ 'excluded'))
)

x2 <- x |> subset(!stratum %in% 'excluded')

x2 <- x2 |> group_by(plot, taxon, stratum) |> summarise(cover = cover.agg(cover), ht.min = mean(ht.min), ht.max = mean(ht.max))
x2 <- x2 |> group_by(plot,taxon) |> mutate(maxcover = max(cover), maxht = max(ht.max))
x2 <- x2 |> subset((ht.max == maxht & cover >= 10) | (maxcover < 10 & maxcover == cover))
x2 <- x2 |> group_by(plot,taxon) |> mutate(maxcover = max(cover), maxht = max(ht.max))
x2 <- x2 |> subset((ht.max == maxht & !taxon %in% "" & !is.na(taxon)))

x2 <- x2 |> group_by(plot,stratum) |> mutate(srank = order(order(-cover,-ht.max)))
x2 <- x2 |> group_by(plot) |> mutate(rank = order(order(-cover,-ht.max)))
x2 <- x2 |> group_by(plot, stratum) |> mutate(maxcover = max(cover), maxht = NULL)
x2 <- x2 |> subset(((maxcover == cover & cover >= 10) | rank <=5) & srank <=3)
x2 <- x2 |> group_by(plot) |> mutate(rank = order(order(-as.numeric(stratum), srank)))

associations <- x2 |> subset(select = c(plot)) |> unique()
plots <- associations$plot
associations$association <- ''
for (i in 1:length(plots)){#i=1
  Com.B <- subset(x2, plot %in% plots[i])
  nrank <- length(unique(Com.B$rank))
  assname <- ""
  for (j in 1:nrank){#j=1
    assname <- ifelse(j == 1,Com.B[Com.B$rank %in% j,]$taxon,
                      ifelse(Com.B[Com.B$rank %in% j,]$stratum == Com.B[Com.B$rank %in% (j-1),]$stratum,
                             paste0(assname, '-',Com.B[Com.B$rank %in% j,]$taxon),paste0(assname, '/',Com.B[Com.B$rank %in% j,]$taxon)))

  }
  associations[associations$plot %in% plots[i],]$association <- assname
}

x3 <- x2 |> left_join(associations)
x3 <- get.structure(veg, simple=FALSE)

#structure

x <- veg |> mutate(stratum = case_when(
    is.na(ht.max) | is.na(type) ~ 'excluded',
    ht.max > 5  & type %in% c('shrub/vine', 'tree')  ~ 'tree',
    type %in% c('tree') ~ 'saplings',
    type %in% c('shrub/vine') ~ 'shrub',
    type %in% c('grass/grasslike') ~ 'grass',
    type %in% c('forb') ~ 'forb',
    !type %in% c('tree', 'shrub/vine', 'forb','grass/grasslike') ~ 'moss',
    TRUE ~ 'excluded'))

  x.ht <- veg |> group_by(plot) |> summarise(ht.max = max(ht.max))

  x2 <- x |> group_by(plot, stratum) |> summarise(cover = cover.agg(cover))
  x3 <- data.frame(plot = unique(x2$plot), tree = 0, sapling = 0, shrub=0, forb=0, grass=0, moss=0)
  x3 <- x3 |> left_join(subset(x2, stratum %in% 'tree')) |> mutate(tree=ifelse(is.na(cover),0,cover), stratum=NULL, cover=NULL)
  x3 <- x3 |> left_join(subset(x2, stratum %in% 'sapling')) |> mutate(sapling=ifelse(is.na(cover),0,cover), stratum=NULL, cover=NULL)
  x3 <- x3 |> left_join(subset(x2, stratum %in% 'shrub')) |> mutate(shrub=ifelse(is.na(cover),0,cover), stratum=NULL, cover=NULL)
  x3 <- x3 |> left_join(subset(x2, stratum %in% 'forb')) |> mutate(forb=ifelse(is.na(cover),0,cover), stratum=NULL, cover=NULL)
  x3 <- x3 |> left_join(subset(x2, stratum %in% 'grass')) |> mutate(grass=ifelse(is.na(cover),0,cover), stratum=NULL, cover=NULL)
  x3 <- x3 |> left_join(subset(x2, stratum %in% 'moss')) |> mutate(moss=ifelse(is.na(cover),0,cover), stratum=NULL, cover=NULL)
  x3 <- x3 |> left_join(x.ht)


  x3 <- x3 |> mutate(structure =
                       case_when(tree+shrub+sapling+forb+grass < 10 ~
                                   case_when(moss < 10 ~ 'barren',
                                             TRUE ~ 'mossland'),
                                 TRUE ~ case_when(tree < 10 ~
                                                    case_when(shrub+sapling < 10 ~
                                                                case_when(grass > forb ~  'grassland',
                                                                          TRUE ~ 'forb meadow'),
                                                              TRUE ~ case_when(tree+shrub+sapling < 75 ~
                                                                                 case_when(sapling > shrub ~ 'woodland regen',
                                                                                           TRUE ~ 'open shrubland'),
                                                                               TRUE ~ case_when(sapling > shrub ~ 'forest regen',
                                                                                                TRUE ~ 'shrub thicket')
                                                              )
                                                    ),
                                                  TRUE ~ case_when(tree < 65 ~
                                                                     case_when(shrub+sapling < 10 ~
                                                                                 case_when(ht.max < 15 ~ 'open low woodland',
                                                                                           ht.max < 30 ~ 'open medium woodland',
                                                                                           TRUE ~ 'open high woodland'),
                                                                               TRUE ~ case_when(sapling > shrub ~
                                                                                                  case_when(shrub+sapling < 75 ~
                                                                                                              case_when(ht.max < 15 ~ 'open regenerating low forest gap',
                                                                                                                        ht.max < 30 ~ 'open regenerating medium forest gap',
                                                                                                                        TRUE ~ 'open regenerating high forest gap'),
                                                                                                            TRUE ~ case_when(ht.max < 15 ~ 'dense regenerating low forest gap',
                                                                                                                             ht.max < 30 ~ 'dense regenerating medium forest gap',
                                                                                                                             TRUE ~ 'dense regenerating high forest gap')),
                                                                                                TRUE ~ case_when(shrub+sapling < 75 ~
                                                                                                                   case_when(ht.max < 15 ~ 'open shrubby low woodland',
                                                                                                                             ht.max < 30 ~ 'open shrubby medium woodland',
                                                                                                                             TRUE ~ 'open shrubby high woodland'),
                                                                                                                 TRUE ~ case_when(ht.max < 15 ~ 'shrubby low woodland',
                                                                                                                                  ht.max < 30 ~ 'shrubby medium woodland',
                                                                                                                                  TRUE ~ 'shrubby high woodland'))
                                                                               )
                                                                     ),
                                                                   TRUE ~ case_when(shrub+sapling < 10 ~
                                                                                      case_when(ht.max < 15 ~ 'open low forest',
                                                                                                ht.max < 30 ~ 'open medium forest',
                                                                                                ht.max < 45 ~ 'open high forest',
                                                                                                ht.max < 60 ~ 'open tall forest',
                                                                                                TRUE  ~ 'open giant forest'),

                                                                                    TRUE ~ case_when(sapling > shrub ~
                                                                                                       case_when(ht.max < 15 ~ 'advance regen low forest',
                                                                                                                 ht.max < 30 ~ 'advance regen medium forest',
                                                                                                                 ht.max < 45 ~ 'advance regen high forest',
                                                                                                                 ht.max < 60 ~ 'advance regen tall forest',
                                                                                                                 TRUE  ~ 'advance regen giant forest'),
                                                                                                     TRUE ~
                                                                                                       case_when(ht.max < 15 ~ 'shrubby low forest',
                                                                                                                 ht.max < 30 ~ 'shrubby medium forest',
                                                                                                                 ht.max < 45 ~ 'shrubby high forest',
                                                                                                                 ht.max < 60 ~ 'shrubby tall forest',
                                                                                                                 TRUE  ~ 'shrubby giant forest')))

                                                  )
                                 )))


  taxa <- veg$taxon
  habit = NA
  gho <- vegnasis::gho
  fill.type <- function(taxa, type=NA){
    x  <-  data.frame(taxa=taxa, type = type) |> mutate(GH0 = '', genus = str_split_fixed(taxa , '[[:blank:]]',3)[,1])
    #first try straight join ----
    x <- x |> left_join(taxon.habits[,c('Scientific.Name','GH')], by = c('taxa'='Scientific.Name'), multiple = 'first')
    x <- x |> mutate(GH0 = ifelse(is.na(GH0)| GH0 %in% "", GH, as.character(GH0)))
    x <- x[,1:4]
    #then try synonym join ----
    x <- x |> left_join(syns[,c('acc','syn')], by=c('taxa'='syn'), multiple = 'first') |> left_join(taxon.habits[,c('Scientific.Name','GH')], by = c('acc'='Scientific.Name'), multiple = 'first')
    x <- x |> mutate(GH0 = ifelse(is.na(GH0)| GH0 %in% "", GH, as.character(GH0)))
    x <- x[,1:4]
    #finally try genus only ----
    x <- x |> left_join(genus.habits, by = c('genus'='genus'), multiple = 'first')
    x <- x |> mutate(GH0 = ifelse(is.na(GH0)| GH0 %in% "", GH, as.character(GH0)))
    x <- x[,1:3]
    x <- x |> mutate(nasis = case_when(
      grepl('^T', GH0) ~ 'tree',
      grepl('^S', GH0)| grepl('^L', GH0) | grepl('^E', GH0) ~ 'shrub/vine',
      grepl('^H.G', GH0)~ 'grass/grasslike',
      grepl('^H', GH0)~ 'forb',
      grepl('^N.B', GH0)~ 'moss',
      grepl('^N.L', GH0)~ 'lichen',
      grepl('^N', GH0)~ 'microbiotic crust',
      TRUE ~ NA
    ))
    nasis <- ifelse(is.na(x$type), x$nasis, x$type)
    return(nasis) }
  #upgrade taxonomy ----
  taxa <- veg$taxon
  harmonize.taxa <- function(taxa){
    x  <-  data.frame(taxa=taxa)
    x <- x |> left_join(syns[,c('acc','syn','ac.binomial')], by=c('taxa'='syn'), multiple = 'first')
    x <- x |> mutate(ac.binomial = ifelse(is.na(ac.binomial), taxa, ac.binomial))
    return(x$ac.binomial)}

  harmonize.taxa(taxa)
  #fill USDA PLANTS Symbols ----
  taxa <- veg$taxon
  fill.usda.symbols <- function(taxa, symbol=NA){
    x  <-  data.frame(taxa=taxa, symbol=symbol)
    x <- x |> left_join(usdaplants[,c('taxon','sym')], by=c('taxa'='taxon'), multiple = 'first')
    x <- x |> mutate(symbol = ifelse(is.na(symbol), sym, symbol))
    return(x$symbol)}

  fill.usda.symbols(taxa)
  #fill nativity ----
  library(vegnasis)
  veg <- clean.veg(vegnasis::nasis.veg)
  taxa <- veg$taxon
  region <- 'Northeast'
  nativity=NA
  fill.nativity <- function(taxa, nativity=NA, region=NA){
  x  <-  data.frame(taxa=taxa, nativity = nativity)
  #first try straight join ----
  x <- x |> left_join(natdat, by = c('taxa'='ac.binomial'), multiple = 'first')
   x <- x |> mutate(nativity0 = case_when(
    region %in% 'Northwest' ~ Northwest,
    region %in% 'Southwest' ~ Southwest,
    region %in% 'NorthCentral' ~ NorthCentral,
    region %in% 'Southcentral' ~ Southcentral,
    region %in% 'Northeast' ~ Northeast,
    region %in% 'Southeast' ~ Southeast,
    region %in% 'Alaska' ~ Alaska,
    region %in% 'Hawaii' ~ Hawaii,
    region %in% 'Caribbean' ~ Caribbean,
    region %in% 'CanadaWest' ~ CanadaWest,
    region %in% 'CanadaEast' ~ CanadaEast,
    region %in% 'Arctic' ~ Arctic,
    region %in% 'Mexico' ~ Mexico,
    TRUE ~ Northwest+Southwest+NorthCentral+Southcentral+Northeast+Southeast+Alaska+Hawaii))
   x <- x |> mutate(nativity = ifelse(is.na(nativity), ifelse(nativity0 > 0,'native','introduced'), nativity))
   x <- x[,1:2]
  #then try synonym join ----
  x <- x |> left_join(syns[,c('acc','ac.binomial','syn')], by=c('taxa'='syn'), multiple = 'first') |> left_join(natdat, by = c('ac.binomial'='ac.binomial'), multiple = 'first')
  x <- x |> mutate(nativity0 = case_when(
    region %in% 'Northwest' ~ Northwest,
    region %in% 'Southwest' ~ Southwest,
    region %in% 'NorthCentral' ~ NorthCentral,
    region %in% 'Southcentral' ~ Southcentral,
    region %in% 'Northeast' ~ Northeast,
    region %in% 'Southeast' ~ Southeast,
    region %in% 'Alaska' ~ Alaska,
    region %in% 'Hawaii' ~ Hawaii,
    region %in% 'Caribbean' ~ Caribbean,
    region %in% 'CanadaWest' ~ CanadaWest,
    region %in% 'CanadaEast' ~ CanadaEast,
    region %in% 'Arctic' ~ Arctic,
    region %in% 'Mexico' ~ Mexico,
    TRUE ~ Northwest+Southwest+NorthCentral+Southcentral+Northeast+Southeast+Alaska+Hawaii))
  x <- x |> mutate(nativity = ifelse(is.na(nativity), ifelse(nativity0 > 0,'native','introduced'), nativity))
  x <- x[,1:2]
    return(x$nativity)}












  remotes::install_github("phytoclast/vegnasis", dependencies = FALSE)
  library(vegnasis)

  veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS=F)

  veg.ca<- veg.raw |> subset(grepl('CA', vegplotid))
  veg.raw <- vegnasis::nasis.veg
  veg <- clean.veg(veg.raw)
  veg$nativity <- NA
  veg <- fill.nativity.df(veg, region='Northwest')
  veg$nativity <- NA
  fill.nativity.df <- function(df, region=NA){
    df$nativity <- fill.nativity(fill.nativity(veg$taxon, region=region, veg$nativity))
    return(df)}



  veg <- veg |> mutate(nativity = fill.nativity(taxa=taxon, region="Southeast", nativity = nativity))

  veg$nativity2 <- fill.nativity(veg$taxon, 'Northwest')
  veg$type <- fill.type(veg$taxon, veg$type)
  saveRDS(veg.ca,'x.RDS')
  vegass <- get.assoc(veg)
  vegstruct.simple <- get.structure(veg)
  vegstruct.complex <- get.structure(veg, simple = FALSE)

  veg$hhh <- get.ht.max(veg$taxon)


  veg <- clean.veg.log(obs, obsspp)

  veg <- veg |> mutate(taxon=harmonize.taxa(veg$taxon, fix=T)) |> fill.type.df() |> fill.nativity.df() |> mutate(symbol = fill.usda.symbols(taxon)) |> fill.hts.df()


  forest <-  veg |> mutate(h = ifelse(diam > 0 & BA > 0, ht.max, NA),
                           d = ifelse(diam > 0 & BA > 0, diam, NA),
                           b = ifelse(diam > 0 & BA > 0, BA, NA)) |>
    group_by(plot) |> filter(ht.max > 5) |>
    summarise(cover = cover.agg(cover), BA = sum(BA, na.rm = T),
              h=sum(h*b, na.rm = T), d=sum(d*b, na.rm = T), b=sum(b, na.rm = T)) |> mutate(BA = ifelse(BA > 0, BA, NA),
                                                                                           ht = ifelse(b == 0, NA, h/b),
                                                                                           diam = ifelse(b == 0, NA, d/b),
                                                                                           h = NULL, d = NULL, b=NULL)
  forest <- forest |> mutate(lcover = log(cover), lht=log(ht), ldiam=log(diam), lBA=log(BA), BA2=BA^2, BA.5=BA^.5,
                             ht2 = ht^2, ht.5 = ht^0.5, diam.5 = diam^0.5, diam2 = diam^2, cover2=cover^2,cover.5=cover^.5)
  library(ggplot2)
  ggplot(forest, aes(x=BA, y=cover))+
    geom_smooth()
  ggplot(forest, aes(x=ldiam, y=lht))+
    geom_smooth()

  mod <- lm(BA ~ cover, forest)
  summary(mod)

  mod <- lm(diam ~ cover+ht+BA, forest)
  step(mod)
  mod <- lm(ldiam ~ lcover+lht+lBA, forest)
  step(mod)

   mod <- lm(diam ~ cover+ht, forest)
   summary(mod)
   mod <- lm(diam ~ ht+ht2+ht.5+lht, forest)
   step(mod)
   mod <- lm(diam ~ lht+ht+ht2, forest)
   step(mod)

   cor(forest[,c('diam','ldiam','diam.5','diam2','ht','lht','ht.5','ht2')], use = 'pairwise.complete.obs')

   ggplot(forest, aes(x=ht.5, y=ldiam))+
     geom_smooth()

   mod <- lm(ldiam ~ ht.5, forest)
   summary(mod)

   h <- c(5, 10, 15, 20, 30, 40)

   d <- exp(0.70799 + 0.61045*h^0.5)
   round(d, 0)

   cor(forest[,c('cover','lcover','cover.5','cover2','BA','lBA','BA.5','BA2')], use = 'pairwise.complete.obs')

   ggplot(forest, aes(x=BA, y=cover))+
     geom_smooth()

   mod <- lm((diam) ~ ht.5, forest)
   summary(mod)

   library(minpack.lm)
   library(growthmodels)
   xy <- forest |> subset(!is.na(BA) & !is.na(cover))
   x <- xy$BA
   y <- xy$cover
   mod1 <- nlsLM(y ~ growthmodels::gompertz(x, alpha, beta, k), start = list(alpha = 100, beta = 1, k = 1))#can swap out start values for fixed


   mod1 <- minpack.lm::nlsLM(y ~ b1*(1-exp(b2*x))^(b3) , start = list(b1=100, b2=-1, b3=1))#can swap out start values for fixed model values
   summary(mod1)

   BA.to.cover <- function(x){
     b1= 94.49218
     b2= -0.08282
     b3= 0.78939
     y = b1*(1-exp(b2*x))^b3
     return(round(y,1))}



   BA=c(0,1,2,5,10,20,50)
  BA.to.cover(BA)

  obssites = vegnasis::obs
  obstaxa = vegnasis::obsspp
  veg <- clean.veg.log(obssites, obstaxa)
  veg <- clean.veg(nasis.veg)
  veg <- veg |> fill.hts.df()

  #Example data created to look as if imported from random csv file.
  obsite <- c('plot1','plot1','plot1', 'plot2', 'plot2')
  obsspp <- c('Acer rubrum','Pinus strobus','Pteridium aquilinum', 'Lindera benzoin', 'Trillium grandiflorum')
  abund <- c(80,10,30,10,10)
  mydata <- data.frame(obsite=obsite, obsspp=obsspp, abund=abund)

  #Identify columns containing data cooresponding to standard column names.
  mydata <- mydata |> mutate(taxon=obsspp, cover=abund, plot=obsite)
  veg <- mydata |> pre.fill.veg()


  veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS = TRUE)
 veg <- clean.veg(veg.raw)

veg <- veg |> fill.nativity.df() |> fill.type.df() |> fill.hts.df()
ass <- get.assoc(veg)
structure <- veg |> get.structure()
structure.alt <- veg |> get.structure.alt()
structure2 <- veg |> get.structure(simple = FALSE)
structure2.alt <- veg |> get.structure.alt(simple = FALSE)
veg <- veg |> summary.strata()

