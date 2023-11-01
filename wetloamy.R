setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
library(sf)
library(terra)
library(soilDB)
library(aqp)

veg.spp <- read.delim('data/Observed_Species.txt', encoding = 'latin1')
veg.site <- read.delim('data/Sites.txt')
veg.site <- subset(veg.site,Latitude != 0 & Observer_Code %in% c('BEL.JH', 'TOL.NB', 'GRR.NJL', 'GRR.GJS') &
                     Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #
veg.site$Observation_Label <- veg.site$Community_Name

veg <- clean.veg.log(veg.site, veg.spp) 

wetloamy <- subset(veg.site, Loamy %in% 'Yes' & Hydric %in% "Yes" )

veg.wetloamy <- subset(veg, plot %in% wetloamy$Observation_ID) |> fill.hts.df() |> fill.type.df() |> fill.nativity.df()

m <- make.plot.matrix(veg.wetloamy, tr = 'log', label='label', taxon='taxon')
d = vegan::vegdist(m, method='bray')
t <- cluster::agnes(d, method = 'ward')|> as.hclust()
k = 7
groups <- cutree(t, k = k)
groups <- dendrogrouporder(t, groups)
a = 'wet loamy vegetation'
folder = 'sitedata'
export.dendro(a,d,t,groups,folder)