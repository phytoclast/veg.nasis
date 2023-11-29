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
veg.site$Observation_Label <- veg.site$Observation_ID

for(i in 1:ncol(veg.site)){
if(!FALSE %in% grepl('Yes|No',veg.site[,i])){
  veg.site[,i] <- ifelse(veg.site[,i] %in% "Yes",1,0)
}}
veg.site <- veg.site |> mutate(Forest = ifelse(Structure %in% 'forest',1,ifelse(Structure %in% 'woodland',0.5,0)))
colnames(veg.site)
unique(veg.site$MLRA)
wetloamy <- subset(veg.site, Loamy %in% 1 & Hydric %in% 1 & Euic %in% 1 & Floodplain %in% 0 & MLRA %in% c('97A','98A1','98A2','99A', '99B'))
# wetloamy <- veg.site[sample(rownames(veg.site),100),]

veg <- clean.veg.log(veg.site, veg.spp)
veg <- subset(veg, cover > 0)
hascover <- veg |> group_by(plot) |> summarize(cover = cover.agg(cover)) |> subset(cover >= 10)
veg <- subset(veg, plot %in% hascover$plot)
veg.wetloamy <- subset(veg, plot %in% wetloamy$Observation_ID) |> fill.hts.df() |> fill.type.df() |> fill.nativity.df()


veg.abiotic <- subset(veg.site, Observation_ID %in% veg.wetloamy$plot,
                      select = c('Observation_ID', "Upper","Middle","Lower","Coastal",
                                           "Floodplain","Inland","Hydric","Nonhydric",
                                           "Aquatic","Wet","Moist","Dry",
                                           "Mucky","Rocky","Sandy","Loamy",
                                           "Calcareous", "Euic", "Dysic", "Salty",
                                           "Fresh", "Natural", "Seminatural", "Cultural",
                                           "Cold", "Cool", "Mild",
                                           "Warm", "Hot", "Humid", "Subhumid",
                                           "Arid", "Microthermal", "Mesothermal", "Megathermal","Forest"))

rownames(veg.abiotic) <- veg.abiotic$Observation_ID
veg.abiotic<- veg.abiotic[,-1]

d.abiotic <- vegan::vegdist(veg.abiotic, method='euclidean')
d.abiotic <- d.abiotic/max(d.abiotic)

m <- make.plot.matrix(veg.wetloamy, tr = 'log', label='plot', taxon='taxon')
d = vegan::vegdist(m, method='bray')
dtab1 <- as.data.frame(as.matrix(d))
dtab2 <- as.data.frame(as.matrix(d.abiotic))
prop=0.5
d2 <- (d*(1-prop)+d.abiotic*(prop))

t <- cluster::agnes(d2, method = 'ward')|> as.hclust()
t$labels
rdf <- data.frame(Observation_ID = t$labels) |> left_join(veg.site[,c("Observation_ID","Community_Name","Site_Type")])
t$labels <- paste(rdf$Community_Name,row(rdf)[,1])

k = 10
groups <- cutree(t, k = k)
groups <- dendrogrouporder(t, groups)
a = 'wet loamy vegetation'
folder = 'sitedata'
export.dendro(a,d,t,groups,folder)
groupdf <- as.data.frame(groups)
rdfgroup <- cbind(rdf,groupdf)
rownames(rdfgroup) = NULL
veg.group <- veg |> left_join(rdfgroup, by=join_by(plot == Observation_ID)) |> subset(!is.na(groups))

rdfgrouped <- rdfgroup |> left_join(wetloamy, by=join_by(Observation_ID == Observation_ID, Community_Name == Community_Name, Site_Type == Site_Type)) |> subset(select = c("Observation_ID", "groups", "Community_Name", "Upper","Middle","Lower","Coastal",
                                                                                                                                                                         "Floodplain","Inland","Hydric","Nonhydric",
                                                                                                                                                                         "Aquatic","Wet","Moist","Dry",
                                                                                                                                                                         "Mucky","Rocky","Sandy","Loamy",
                                                                                                                                                                         "Calcareous", "Euic", "Dysic", "Salty",
                                                                                                                                                                         "Fresh", "Natural", "Seminatural", "Cultural",
                                                                                                                                                                         "Cold", "Cool", "Mild",
                                                                                                                                                                         "Warm", "Hot", "Humid", "Subhumid",
                                                                                                                                                                         "Arid", "Microthermal", "Mesothermal", "Megathermal","Forest"))
rdfgrouped <- rdfgrouped |> group_by(groups) |> summarise(across(c(Upper,Middle,Lower,Coastal,
                                                Floodplain,Inland,Hydric,Nonhydric,
                                                Aquatic,Wet,Moist,Dry,
                                                Mucky,Rocky,Sandy,Loamy,
                                                Calcareous, Euic, Dysic, Salty,
                                                Fresh, Natural, Seminatural, Cultural,
                                                Cold, Cool, Mild,
                                                Warm, Hot, Humid, Subhumid,
                                                Arid, Microthermal, Mesothermal, Megathermal,Forest), mean))


#grouped summary
#
veg.summary <- veg.group |> fill.nativity.df() |> fill.type.df() |> fill.hts.df() |> mutate(symbol = fill.usda.symbols(taxon))

veg.summary2  <- veg.summary |>  summary.ESIS(group='groups', lowerQ = 0.5, upperQ = 0.95, normalize = TRUE, breaks = c(0.5, 2, 5, 15))


coverbounds <- function(cover.Low, cover.High, cover.min, cover.max, cover.mean){
  upper.x = cover.High + (cover.mean - (cover.High+cover.Low)/2)
  lower.x = cover.Low + (cover.mean - (cover.High+cover.Low)/2)
  upper.x = upper.x + (lower.x - pmax(lower.x, cover.min))
  lower.x = pmax(lower.x, cover.min)
  lower.x = lower.x + (upper.x - pmin(upper.x, cover.max))
  upper.x = pmin(upper.x, cover.max)
  return(c(lower.x,upper.x))}

#example numbers
cover.min = 10
cover.Low = 27.5
cover.mean = 33.8
cover.High = 74.0
cover.max = 80.0
x = coverbounds(cover.Low, cover.High, cover.min, cover.max, cover.mean)
x
mean(x)

#species with frq of 0.444
cover.min = 0
cover.50 = 0.0
cover.mean = 9
cover.95 = 10
cover.max = 12
x = coverbounds(cover.50, cover.95, cover.min, cover.max, cover.mean)
x
mean(x)
