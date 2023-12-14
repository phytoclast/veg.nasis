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
veg.soil <- subset(veg.site, select = c(Observation_ID, State, County, FIPS, Site_Type, Soil.Taxon, Soil.Series, Map.Unit))

veg <- clean.veg.log(veg.site, veg.spp)
veg <- subset(veg, cover > 0)

# x=veg |> fill.type.df() |> fill.hts.df();simple=TRUE
hascover <- veg |> group_by(plot) |> summarize(cover = cover.agg(cover)) |> subset(cover >= 10)
veg <- subset(veg, plot %in% hascover$plot)
veg <- subset(veg, plot %in% wetloamy$Observation_ID) |> fill.hts.df() |> fill.type.df() |> fill.nativity.df() |>
  mutate(taxon = harmonize.taxa(taxon, fix=T, sensu = 'bonap'),
         taxon = extractTaxon(taxon, 'binomial'), #use for analysis
         usda = harmonize.taxa(taxon, fix=T, sensu = 'usda'),
         symbol = fill.usda.symbols(usda)) #for entry into EDIT
veg.structure <- veg |> get.structure(simple=FALSE)
veg.association <- veg |> get.assoc() |> left_join(veg.structure)
veg.association <- veg.association |> left_join(data.frame(plot = veg.site$Observation_ID, wet = veg.site$Wet, Natural = veg.site$Natural))

veg.association <- veg.association |> 
  mutate(groups= case_when(tree >= 25 & wet %in% 1 & Natural %in% 1 ~ '1.1 Swamp Forest',
                           tree >= 25 & wet %in% 0 & Natural %in% 1 ~ '5.1 Native Drained Forest',
                           tree >= 25 & wet %in% 1 & Natural %in% 0 ~ '3.2 Exotic Drained Forest',
                           tree >= 25 & wet %in% 0 & Natural %in% 0 ~ '4.2 Exotic Swamp Forest',
                           (grepl('shrubland', structure)|grepl('thick', structure)) & wet %in% 1 & Natural %in% 1 ~ '1.5 Inundated Shrub Swamp',
                           tree < 25 & wet %in% 1 & Natural %in% 1 ~ '1.2 Wet Meadow',
                           tree < 25 & wet %in% 0 & Natural %in% 0 ~ '3.1 Exotic Drained Meadow & Shrub',
                           tree < 25 & wet %in% 1 & Natural %in% 0 ~ '4.1 Exotic Wet Meadow & Shrub',
                           tree < 25 & wet %in% 0 & Natural %in% 1 ~ '5.2 Native Drained Meadow & Shrub',
                           TRUE ~ 'other'))




veg.abiotic <- subset(veg.site, Observation_ID %in% veg$plot)
veg.abiotic<- veg.abiotic |> left_join(veg.structure, by=join_by(Observation_ID==plot)) |> 
  subset(select = c('Observation_ID', "Wet","Moist","Natural", "Seminatural", "Forest")) 
# subset(select = c('Observation_ID', "Wet","Moist","Natural", "Seminatural", "tree","sapling","shrub","forb","grass","ht.max")) 
# select = c('Observation_ID', "Upper","Middle","Lower","Coastal",
#            "Floodplain","Inland","Hydric","Nonhydric",
#            "Aquatic","Wet","Moist","Dry",
#            "Mucky","Rocky","Sandy","Loamy",
#            "Calcareous", "Euic", "Dysic", "Salty",
#            "Fresh", "Natural", "Seminatural", "Cultural",
#            "Cold", "Cool", "Mild",
#            "Warm", "Hot", "Humid", "Subhumid",
#            "Arid", "Microthermal", "Mesothermal", "Megathermal","Forest"))

rownames(veg.abiotic) <- veg.abiotic$Observation_ID
veg.abiotic<- veg.abiotic[,-1]

d.abiotic <- vegan::vegdist(veg.abiotic, method='euclidean')
d.abiotic <- d.abiotic/max(d.abiotic)

m <- make.plot.matrix(veg, tr = 'log', label='plot', taxon='taxon')
d = vegan::vegdist(m, method='bray')
dtab1 <- as.data.frame(as.matrix(d))
dtab2 <- as.data.frame(as.matrix(d.abiotic))
prop=0.9
d2 <- (d*(1-prop)+d.abiotic*(prop))

t <- cluster::agnes(d2, method = 'ward')|> as.hclust()
t$labels
rdf <- data.frame(Observation_ID = t$labels) |> left_join(veg.site[,c("Observation_ID","Community_Name","Site_Type")])
t$labels <- paste(rdf$Community_Name,row(rdf)[,1])

k = 11-3
groups <- cutree(t, k = k)
groups <- dendrogrouporder(t, groups)
a = 'wet loamy vegetation'
folder = 'sitedata'
export.dendro(a,d,t,groups,folder)
groupdf <- as.data.frame(groups)
rdfgroup <- cbind(rdf,groupdf)
rownames(rdfgroup) = NULL

#link cluster groups
veg.group <- veg |> left_join(rdfgroup, by=join_by(plot == Observation_ID)) |> subset(!is.na(groups))


#link manually determined groups based on veg structure
veg.group <- veg |> left_join(veg.association[,c('plot','groups')], by=join_by(plot == plot)) |> subset(!is.na(groups))



rdfgrouped <- rdfgroup |> left_join(wetloamy, by=join_by(Observation_ID == Observation_ID, Community_Name == Community_Name, Site_Type == Site_Type)) |> subset(select = c("Observation_ID", "groups", "Community_Name", "Upper","Middle","Lower","Coastal",                                                                     "Floodplain","Inland","Hydric","Nonhydric",
                                                                                                                                                                           "Aquatic","Wet","Moist","Dry",
                                                                                                                                                                           "Mucky","Rocky","Sandy","Loamy",                                                                                                             "Calcareous", "Euic", "Dysic", "Salty",                                                                                                                    "Fresh", "Natural", "Seminatural", "Cultural",                                                                                                              "Cold", "Cool", "Mild",
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
mlra99plots <- subset(veg.site, MLRA %in% c('99A', '99B'))
mlra97plots <- subset(veg.site, MLRA %in% c('97A'))
mlra98A2plots <- subset(veg.site, MLRA %in% c('98A2'))

theseplots = mlra99plots$Observation_ID #plots to emphasize
thoseplots = c(mlra97plots$Observation_ID,mlra98A2plots$Observation_ID) #plots to de-emphasize
# theseplots = mlra97plots$Observation_ID
# thoseplots = mlra99plots$Observation_ID
groupcounts <-  veg.group |> mutate(these = ifelse(plot %in% theseplots, 1,0), those = ifelse(plot %in% thoseplots, 1,0)) |> subset(select=c(groups, plot, these, those)) |> unique() |> group_by(groups) |> summarise(nplot = length(plot), these = sum(these), those = sum(those))
groupcounts <- groupcounts |> mutate(prewt1 = 1/(nplot+1), prewt2 = 11/(nplot-those+1),prewt3 = 11/(these+1),actual = these/nplot,
                                     strength = prewt3*these/(prewt3*these + prewt1*those + prewt2*(nplot-these-those)))
veg.summary <- veg.group  |> left_join(groupcounts)|> mutate(wt = ifelse(plot %in% theseplots, prewt3,prewt2),
                                                             wt = ifelse(plot %in% thoseplots,prewt1,wt))

subset(veg.group, plot %in% mlra99plots$Observation_ID, select=c(groups, plot)) |> unique()



EDIT  <- veg.summary |>  summary.ESIS.wt(group='groups', wt='wt', lowerQ = 0.5, upperQ = 0.95, normalize = TRUE, breaks = c(0.5, 2, 5, 15), forEDIT = F)



group_community <- data.frame(plot = EDIT$group,
                              taxon = EDIT$taxon,
                              symbol = EDIT$symbol,
                              type = EDIT$type,
                              nativity=EDIT$nativity,
                              cover=EDIT$cover.mean,
                              stratum.min = EDIT$stratum.min,
                              stratum.max = EDIT$stratum.max,
                              ht.min = EDIT$Bottom,
                              ht.max = EDIT$Top) |> pre.fill.veg()


group_community.ass <- group_community |> fill.hts.df()
group_community.ass <- get.assoc(group_community.ass)
# group_community.ass <- group_community.ass |> left_join(rdfgrouped[,c('groups','Wet','Forest','Natural')], by = join_by(plot==groups))


#structural summary ---- 
unique(veg.group$groups)
veg.str <- subset(veg.group, groups %in% "4.2 Exotic Swamp Forest") |> summary.crown.thickness(breaks = c(0.15,0.3,0.6,1.4,4,12,24,37)) |> structure.fill.zero() |> 
                                                                subset(type %in% c('tree', 'shrub/vine', 'grass/grasslike',  'forb'))
veg.str.pct <- veg.str |> group_by(type, stratum, stratum.label, bottom, top) |>
  summarise(X25 = quantile(Cover, 0.25),
            X75 = quantile(Cover, 0.75))
tree <- subset(veg.str.pct, type %in% 'tree', select=c(stratum, X25, X75)) ; colnames(tree) <- c('stratum','t25','t75')
shrub <- subset(veg.str.pct, type %in% 'shrub/vine', select=c(stratum, X25, X75)) ; colnames(shrub) <- c('stratum','s25','s75')
grass <- subset(veg.str.pct, type %in% 'grass/grasslike', select=c(stratum, X25, X75)) ; colnames(grass) <- c('stratum','g25','g75')
forb <- subset(veg.str.pct, type %in% 'forb', select=c(stratum, X25, X75)) ; colnames(forb) <- c('stratum','f25','f75')

veg.str.pct1 <- unique(veg.str.pct[,2:3]) |> left_join(tree) |> left_join(shrub) |> left_join(grass) |> left_join(forb)

#composition summary ----
EDIT  <- veg.summary |>  summary.ESIS.wt(group='groups', wt='wt', lowerQ = 0.5, upperQ = 0.95, normalize = TRUE, breaks = c(0.5, 5, 15), forEDIT = T)

EDIT <- subset(EDIT, cover.High >=0.1 & frq.plot >= 0.05)
write.csv(EDIT, 'wetloamy.csv', row.names = F, na="")
