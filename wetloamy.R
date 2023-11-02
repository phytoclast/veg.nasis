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
# wetloamy <- subset(veg.site, Loamy %in% 1 & Hydric %in% 1 & Floodplain %in% 0)
wetloamy <- veg.site[sample(rownames(veg.site),100),]
veg <- clean.veg.log(veg.site, veg.spp) 

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
t$labels <- paste(rdf$Site_Type,row(rdf)[,1])

k = 7
groups <- cutree(t, k = k)
groups <- dendrogrouporder(t, groups)
a = 'wet loamy vegetation'
folder = 'sitedata'
export.dendro(a,d,t,groups,folder)