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
# nasis.veg <- readRDS('data_raw/veg.raw.select.RDS')
# usethis::use_data(nasis.veg, overwrite = T)



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

veg.raw <-  vegnasis::nasis.veg
veg <- clean.veg(veg.raw)

vegs <- veg |> subset(plot %in% '2022MI165023.P') |> fill.hts.df()

vegs <- vegs |> mutate(diam =  fill.diameters(ht.max,diam))
vegs <- vegs |> mutate(stems = trees_per_ha(BA, diam))

totalBA <- sum(vegs$BA, na.rm = T)
vegs.sum <- vegs |> mutate(cancov = ifelse(ht.max > 5, cover, NA)) |> group_by(plot) |> mutate(cancov.sum = sum(cancov, na.rm = T),
                                                                                               BA.sum=sum(BA, na.rm = T),
                                                                                               BA2 = BA.sum*cancov/cancov.sum)

