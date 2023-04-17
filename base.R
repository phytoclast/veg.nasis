library(soilDB)
library(stringr)
library(dplyr)

#install and load package
remotes::install_github("phytoclast/vegnasis", dependencies = FALSE)
library(vegnasis)
#fresh NASIS data
veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS = F)
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
# veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS = F)
# select1 <- subset(nasis.veg, !grepl('2022MI', vegplotid))
# select2 <- subset(veg.raw, grepl('2022MI', vegplotid))
# select3 <- rbind(select2,select1)
# saveRDS(select3,'data_raw/veg.raw.select.RDS')

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
library(vegnasis)
veg.raw <-  vegnasis::nasis.veg
veg <- clean.veg(veg.raw) |> fill.type.df() |> fill.hts.df()
veg$habit <- get.habit.code(veg$taxon)
veg.s <- subset(veg,  grepl('2022MI165021.P',plot))
taxon <- c('Acer rubrum', 'Pinus resinosa')
crfill <- c(NA,"#80991A")
stfill <- c('gray',"#B36666")
crshape <- c(NA,'conifer2')
override <- data.frame(taxon=taxon,stfill=stfill,crfill=crfill,crshape=crshape)
veg.s <- veg.s |> left_join(override)

plants <- grow_plants(veg.s, pwidth=50)

veg_profile_plot0(plants, xslope=50, yslope=25)
'2022MI165021.P'
'2021WA031024'
rgb(0.7,0.4,0.4)
#Create fake data
plot <- c('plot1')
taxon <- c('Quercus macrocarpa','UNK', 'Festuca', 'Andropogon', 'Liatris')
cover <- c(20,20,60,10,5)
crown.max <- c(15, 30, 0.6, 2, 10)
dbh <- c(60,NA,NA,NA,NA)
habit <- c(NA,'T.NE',NA, NA,'S.BD')
crfill <- c('yellow',NA,'red',NA,NA)
crcolor <- c('green','purple','red',NA,NA)
stfill <- c('white',NA,'white',NA,NA)
stcolor <- c('gray','blue','red',NA,NA)

mydata <- data.frame(taxon=taxon, cover=cover, plot=plot, crown.max = crown.max, dbh.max = dbh,
                     habit=habit, crfill=crfill, crcolor=crcolor, stfill=stfill, stcolor=stcolor)

veg <- mydata |> pre.fill.veg()
plants <- grow_plants(veg)
veg_profile_plot(plants, unit='m',  skycolor = rgb(0.8,0.98,1), fadecolor = 'lightgray', gridalpha = 0.0, groundcolor = rgb(0.55,0.45,0.2), ylim = c(0,20))


veg.raw <- readRDS('data/nasisvegraw3.RDS')
veg <- clean.veg(veg.raw) |> fill.hts.df()
ass <- veg |> get.assoc()
veg <- veg |> left_join(ass)

m <- veg |> make.plot.matrix(tr = 'log', rc = T, nr=T, label='association')
#distance matrix based on Bray-Curtis simularity.
d = vegan::vegdist(m, method='bray')
#Cluster analysis using Ward's method using distance matrix.
t <- cluster::agnes(d, method = 'ward')|> as.hclust()
#Define number of groups to color the dendrogram by.
k = 4
groups <- cutree(t, k = k)
#This function rearranges the branchs and groups so that the tree is always oriented with most nested branches to the bottom of the plot (when tree oriented vertically with branches to the right).
groups <- dendrogrouporder(t, groups)
a = 'Vegetation of Michigan'
plot.dendro(a,d,t,groups)



hb <- vegnasis::taxon.habits
s <- vegnasis::shapes




plot <- c('plot1')
taxon <- c('Quercus macrocarpa','UNK','Pteridium', 'Festuca', 'Andropogon', 'Liatris')
cover <- c(20,5,10,60,10,5)
crown.max <- c(15,4,1,0.6,2,0.4)
crfill <- c(NA,"#99E6B3",NA,NA,NA,NA)
dbh <- c(60,NA,NA,NA,NA,NA)
habit <- c(NA,'S.BD',NA,NA,NA,NA)
mydata <- data.frame(plot=plot,taxon=taxon, cover=cover, habit=habit, crown.max = crown.max, dbh.max = dbh,crfill=crfill)


veg <- mydata |> pre.fill.veg()
plants <- grow_plants(veg, plength=100) #Grow more plants on a longer 100 m plot (default was 50 m).
veg_profile_plot(plants, unit='m',  skycolor = rgb(0.8,0.98,1), fadecolor = 'lightgray', gridalpha = 0.1, groundcolor = rgb(0.55,0.45,0.2), xlim=c(0,100))

rgb(0.6,0.9,0.7)






veg <- veg  |> fill.hts.df()|> subset(grepl('WA',plot))

veg <- veg |> mutate(dbh.r =  fill.diameters(ht.max,dbh.max,dbh.min))
veg <- veg |> mutate(cw =  case_when(type %in% 'tree' | ht.max > 5 ~ pmax(est_crown_width(dbh.r),1),
                                     type %in% 'shrub/vine' ~ pmax(pmin(3,ht.max),1),
                                     TRUE ~ 1))
veg <- veg |> mutate(density =  density_from_cw(cover, cw))
veg <- veg |> mutate(BA.r =  BA_per_ha(density, dbh.r))

veg <- veg |> group_by(plot) |> mutate(BA.sum = sum(BA, na.rm = T), BA.rsum = sum(BA.r, na.rm = T), BA.sum = ifelse(is.na(BA.sum), BA.rsum,BA.sum), BA.ratio = BA.sum/BA.rsum,  BA.rsum = NULL)#

veg <- veg |> mutate(BA.r = ifelse(ht.max <= 5, BA.r, round(BA.r*BA.ratio,1)),
                     density = ifelse(ht.max <= 5, density, round(density*BA.ratio,0)),
                     cw = ifelse(ht.max <= 5, cw, round(cw*BA.ratio^-0.5,1)))

veg <- veg |> mutate(habit= get.habit.code(taxon),
                     crshape = case_when(grepl('^T', habit) & grepl('N', habit) ~ 'conifer1',
                                         grepl('^T', habit)  ~ 'blob',
                                         type %in% 'shrub/vine' ~ 'cloud1',
                                         grepl('FE', habit) ~ 'ferny',
                                         grepl('F', habit) ~ 'forby',
                                         type %in% 'grass/grasslike' ~ 'grassy'),
                     stshape = case_when(grepl('^T', habit) & grepl('N', habit) ~ 'trunk',
                                         grepl('^T', habit)  ~ 'trunk',
                                         type %in% 'shrub/vine' ~ 'sticks',
                                         grepl('FE', habit) ~ NA,
                                         grepl('F', habit) ~ NA,
                                         type %in% 'grass/grasslike' ~ NA),
                     fun=case_when(grepl('^T', habit)  ~ 'T',
                                   grepl('^S', habit)  ~ 'S',
                                   grepl('^H', habit)  ~ 'H'),
                     stems = round(0.1*density,0) #count number of stem objects required for tenth hectare plot
                     )
strats <- veg |> subset(fun %in% c('T','S','H') & stems > 0) |> arrange(plot,-ht.max, -cover)
strats$seq <- c(1:nrow(strats))
strats <- strats |> group_by(plot) |> mutate(seqmin = min(seq), seq = seq-seqmin+1, seqmin = NULL)


strats <- subset(strats,  grepl('2021WA031024',plot))
#make stand
stand <- make_hex_stand(0.5,1) |> subset(yp >= 15 & yp < 35) |> mutate(wtn = wt, stratid = NA)
#define counts per stratum





for (i in 1:nrow(strats)){#i=1
  thistrat = strats$seq[i]
  nstems = strats$stems[i]
  newstumps <- sample(stand$stumpid, size = nstems, prob = stand$wtn, replace = T)
  stand <- stand |> mutate(wtn = ifelse(stand$stumpid %in% newstumps, 0, wtn),
                           stratid = ifelse(stand$stumpid %in% newstumps, thistrat,stratid))
}

#Create shapes of the right size, then distribute into the stump positions.
for (i in 1:nrow(strats)){#i=1
  thistrat <- strats[i,]
  plant0 <- make_plant(thistrat$fun, thistrat$ht.max, thistrat$ht.min,thistrat$cw,thistrat$dbh.r, thistrat$crshape, thistrat$stshape)
  stumps0 <- stand |> subset(stratid %in% i)
  plant0 <- merge(stumps0, plant0) |> mutate(objid = paste0(stratid,obj,stumpid))
  if(i==1){plants <- plant0}else{plants <- rbind(plants,plant0)}
}

#randomize sizes and positions
plants <- plants |> group_by(stumpid) |>
  mutate(ht.max = max(z), crwd = max(x)-min(x),
         xpp = xp + runif(1, min = -0.8, max = 0.8),#shift position on grid
         zr = rnorm(1,ht.max, ht.max/10)/ht.max,#deviation in height
         xr = (rnorm(1,ht.max, ht.max/10)/ht.max+rnorm(1,crwd, crwd/10)/crwd)/2,#deviation in width partially related to height
         xn = x*xr+xpp,#resized width and put on new position
         zn = z*zr*(1-1/15))#resized height adjusted downward show that variation is less than max height in the field

#rearrange stems depth drawing order
plants <- plants |> arrange(yp,stumpid, objid, ptord)
ypmax <- max(plants$yp)
ypmin <- min(plants$yp)
ypwid <- ypmax-ypmin
plants <- plants |> mutate(depth = case_when(yp < ypmin+ypwid*(1/5) ~ 'E',
                                             yp < ypmin+ypwid*(2/5) ~ 'D',
                                             yp < ypmin+ypwid*(3/5) ~ 'C',
                                             yp < ypmin+ypwid*(4/5) ~ 'B',
                                             TRUE ~ 'A'))

crowns1 <- plants |> subset(depth %in% 'E' & obj %in% c('crown','herb')) |>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.8), color=colormixer(color, "#D9F2FF", 0.8))
crowns2 <- plants |> subset(depth %in% 'D' & obj %in% c('crown','herb')) |>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.6), color=colormixer(color, "#D9F2FF", 0.6))
crowns3 <- plants |> subset(depth %in% 'C' & obj %in% c('crown','herb')) |>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.4), color=colormixer(color, "#D9F2FF", 0.4))
crowns4 <- plants |> subset(depth %in% 'B' & obj %in% c('crown','herb'))|>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.2), color=colormixer(color, "#D9F2FF", 0.2))
crowns5 <- plants |> subset(depth %in% 'A' & obj %in% c('crown','herb'))

stems1 <- plants |> subset(depth %in% 'E' & obj %in% 'stem')|>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.8), color=colormixer(color, "#D9F2FF", 0.8))
stems2 <- plants |> subset(depth %in% 'D' & obj %in% 'stem')|>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.6), color=colormixer(color, "#D9F2FF", 0.6))
stems3 <- plants |> subset(depth %in% 'C' & obj %in% 'stem')|>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.4), color=colormixer(color, "#D9F2FF", 0.4))
stems4 <- plants |> subset(depth %in% 'B' & obj %in% 'stem')|>
  mutate(fill=colormixer(fill, "#D9F2FF", 0.2), color=colormixer(color, "#D9F2FF", 0.2))
stems5 <- plants |> subset(depth %in% 'A' & obj %in% 'stem')

plants2 <- rbind(crowns1,crowns2,crowns3,crowns4,crowns5,stems1,stems2,stems3,stems4,stems5)

# ggplot()+
#   geom_point(aes(x=1,y=1))+
#   theme(panel.background = element_rect(fill = rgb(0.4,0.6,0.5)))

pcolor <- plants2$color |> unique() |> sort()
pfill <- plants2$fill |> unique()|> sort()



ggplot() +
  geom_polygon(data=stems1, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=crowns1, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=stems2, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=crowns2, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=stems3, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=crowns3, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=stems4, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=crowns4, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=stems5, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  geom_polygon(data=crowns5, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
  scale_fill_manual(values=pfill)+
  scale_color_manual(values=pcolor)+
  theme(legend.position = "none",
        panel.background = element_rect(fill = rgb(0.85,0.95,1,0.5),
                                        colour = "black",
                                        linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                        colour = rgb(0.1, 0.1, 0.1, 0.3)),
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid',
                                        colour = rgb(0.1, 0.1, 0.1, 0.1))
  )+
  coord_fixed(ratio = 1)+
  scale_y_continuous(trans='identity', breaks = c(-10:(120/5))*5,minor_breaks = c(-10:(120)), limits = c(0,100))+
  scale_x_continuous(breaks = c(-10:(120/5))*5, minor_breaks = c(-10:(120)), limits = c(-10,60))


































vegs <- vegs |> mutate(stratum = case_when(ht.max > 45 ~ 'T3',
                                           ht.max > 15 ~ 'T2',
                                           ht.max > 5 & type %in% c('tree', 'shrub/vine') ~ 'T1',
                                           ht.max > 0.5 & type %in% c('tree', 'shrub/vine') ~ 'S',
                                           ht.max >= 0 & !type %in% c('moss', 'lichen') ~ 'H',
                                           type %in% c('moss', 'lichen') ~ 'M',
                                           ht.max < 0 ~ 'A',
                                           ))

vegs <- vegs |> group_by(plot, stratum) |> mutate(stratum.cover = cover.agg(cover))


veg <- clean.veg.log(obs, obsspp)

veg <- veg |> mutate(taxon=harmonize.taxa(veg$taxon, fix=T)) |> fill.type.df() |> fill.nativity.df() |> mutate(symbol = fill.usda.symbols(taxon)) |> fill.hts.df()
forest <-  veg |> mutate(diam = ifelse(is.na(dbh.min),dbh.max,(dbh.max+dbh.min)/2),
                         h = ifelse(diam > 0 & BA > 0, ht.max, NA),
                         d = ifelse(diam > 0 & BA > 0, diam, NA),
                         b = ifelse(diam > 0 & BA > 0, BA, NA)) |>
  group_by(plot) |> filter(ht.max > 5) |>
  summarise(cover = cover.agg(cover), BA = sum(BA, na.rm = T),
            h=sum(h*b, na.rm = T), d=sum(d*b, na.rm = T), b=sum(b, na.rm = T)) |> mutate(BA = ifelse(BA > 0, BA, NA),
                                                                                         ht = ifelse(b == 0, NA, h/b),
                                                                                         diam = ifelse(b == 0, NA, d/b),
                                                                                         h = NULL, d = NULL, b=NULL)
n <- c(0:100)*400
(1/10000)*10000*100
c0 <- (1/10000)*n*100
c1 <- pmin((1/10000)*n*100, 100)
c2 <- (1-(1-1/10000)^n)*100
cf = ((100-c2)/1+0)/100+0.2
c3 <- c1*2/3+c2*1/3

w = 10
c = 90
cf=0.8
n1 <- log(1-c/100)/log(1-w/10000)
n2 <- c/(100*(w/10000))
n3 <- n1*(1-cf)+n2*cf
c1 <- pmin((w/10000)*n3*100, 100)
c2 <- (1-(1-w/10000)^n3)*100
c1*2/3+c2*1/3

ht <- c(0.5,2,5,10,20,30,45,80)
dbh <- fill.diameters(ht)
dbh
dbh <-exp(1.248705*log(ht)) |> round(1)
dbh
cover <- c(10)
BA = cover/100*30
density = trees_per_ha(BA, dbh)
cw = est_crown_width(density, cover, dbh)
cw
#45.2 m tall, 134.7 cm


# treedb <- read.csv('C:/workspace2/vegrob/All Trees (m).csv')
#
# treedb <- treedb |> mutate(ht.max =  readr::parse_number(Height..m.),
#                            dbh =  readr::parse_number(Girth..cm.)/3.141592,
#                            cw = readr::parse_number(Crown.spread..m.)) |> subset(select = c(Botanical.Name, ht.max,dbh,cw))
# write.csv(treedb, 'data/treedb.csv', row.names = F)
treedb <- read.csv('data/treedb.csv') |> subset(dbh < 900 & ht.max > 5 & dbh > 10)

forest <- treedb |> mutate(ht=ht.max, lht=log(ht), ldiam=log(dbh),
                           ht2 = ht^2, ht.5 = ht^0.5, diam.5 = dbh^0.5, diam2 = dbh^2, lcw = log(cw)) |> subset(!is.na(ht) & !is.na(dbh))
#To test bends in the curve, several transformations to see which fit was best.
cor(forest[,c('dbh','ldiam','diam.5','diam2','ht','lht','ht.5','ht2','cw','lcw')], use = 'pairwise.complete.obs')

mod <- lm(log(dbh)~log(ht.max)+0, data=forest)
summary(mod)
library(ggplot2)
forest <- forest |> mutate(dbhmod = exp(1.248705*log(ht.max)), dbhold = fill.diameters(ht.max))
ggplot(data=forest)+
  geom_point(aes(x=ht.max, y=dbh), color='black', size= 0.1, alpha=0.1)+
  geom_smooth(aes(x=ht.max, y=dbh), color='red')+
  geom_smooth(aes(x=ht.max, y=dbhold), color='blue')+
  geom_smooth(aes(x=ht.max, y=dbhmod), color='green')+
  scale_x_continuous(breaks = c(0:20)*20)+#, limits = c(0,200))+
  scale_y_continuous(breaks = c(0:40)*50)#, limits = c(0,100))

mod <- lm(cw~dbh+I(dbh^2), data=forest)
summary(mod)

mod1 <- minpack.lm::nlsLM(cw ~ b1*(1-exp(b2*dbh))^(b3), data=forest , start = list(b1=100, b2=-1, b3=1))
summary(mod1)
b1= 69.690000
b2= -0.002392
b3=  0.857591
forest <- forest |> mutate(cwmod = b1*(1-exp(b2*dbh))^(b3), ratio=cw/dbh, mratio=cwmod/dbh)
ggplot(data=forest)+
  geom_point(aes(x=dbh, y=cw), color='black', size= 0.1, alpha=0.5)+
  geom_smooth(aes(x=dbh, y=cw), color='red')+
  geom_smooth(aes(x=dbh, y=cwmod), color='blue')+

  scale_x_continuous(breaks = c(0:20)*50)+#, limits = c(0,200))+
  scale_y_continuous(breaks = c(0:40)*20)#, limits = c(0,100))
ggplot(data=forest)+
  geom_point(aes(x=dbh, y=ratio), color='black', size= 0.1, alpha=0.5)+
  geom_smooth(aes(x=dbh, y=ratio), color='red')+
  geom_smooth(aes(x=dbh, y=mratio), color='blue')+

  scale_x_continuous(breaks = c(0:20)*50)+#, limits = c(0,200))+
  scale_y_continuous(breaks = c(0:40)*.20)#, limits = c(0,100))



cover = 90
cw = 1
density_from_cw(cover, cw)


c1 <- pmin((w/10000)*n2*100, 100)
c2 <- (1-(1-w/10000)^n2)*100
c1*2/3+c2*1/3

df <-  data.frame(n=n,c0=c0,c1=c1,c2=c2,c3=c3)

library(ggplot2)

ggplot(data=df)+
  geom_line(aes(x=c0, y=c1), color='red')+
  geom_line(aes(x=c0, y=c2), color='blue')+
  geom_line(aes(x=c0, y=c3), color='green')+
  scale_x_continuous(breaks = c(0:20)*20, limits = c(0,200))+
  scale_y_continuous(breaks = c(0:40)*20, limits = c(0,100))



BA = forest$BA
cc = forest$cover

BA=c(0:600)/5
cc=BA.to.cover(BA)
BA2=30*cc/100
df = data.frame(BA=BA,
                cc=cc,
                BA2=BA2)
library(ggplot2)

ggplot(data=forest, aes(x=cover, y=BA))+
  geom_point()+
  geom_smooth()+
  geom_line(data=df, aes(x=cc, y=BA), color='red')+
  geom_smooth(data=df, aes(x=cc, y=BA2), color='green')+
  scale_x_continuous(breaks = c(0:20)*20)+
  scale_y_continuous(breaks = c(0:40)*20)













shape = branch
bnarrow <- function(shape){
  df=NULL
  for(i in 1:50){#i=10

    yu = (i+3)/50
    yl = (i-3)/50
    ys <- subset(shape, y>= yl & y<= yu)
    if(nrow(ys)>0){
      w <- max(ys$x)-min(ys$x)
      y <- (yu+yl)/2

      df0 <- data.frame(w=w,y=y)
      if(is.null(df)){df = df0}else{df=rbind(df,df0)}
    }}

  wmax <- max(df$w)
  ywidest <- mean(subset(df, w %in% wmax)$y)
  wmin <- min(subset(df, y < ywidest)$w)
  ynarrowest <-  mean(subset(df, w %in% wmin)$y)

}


#convert raster image to vector
library(terra)
library(sf)
maple <- rast('C:/workspace2/vegrob/maple.png')
maple <- ifel(maple > 0, NA,1)
maple1 <- as.polygons(maple) |> geom() |> as.data.frame()
ggplot(maple1, aes(x,y,group=paste(part,hole), fill=hole))+
  geom_polygon()
maple1 <- subset(maple1, hole <=0)
maple1 <- maple1 |> group_by(part) |> mutate(size = length(part))
maxpart <- max(maple1$size)
maple1 <- subset(maple1, size %in% maxpart)
ggplot(maple1, aes(x,y,group=paste(part,hole), fill=paste(part,hole)))+
  geom_polygon()




#morph shape top and bottom
shape = maple1
shape$x <- (shape$x - (max(shape$x) + min(shape$x))/2) / (max(shape$x) - min(shape$x))
shape$y <- 1-1*(shape$y - min(shape$y)) / (max(shape$y) - min(shape$y))
shape$y <- shape$y*-1+1
df=NULL
for(i in 1:25){#i=10

  yu = (i+3)/25
  yl = (i-3)/25
  ys <- subset(shape, y>= yl & y<= yu)
  if(nrow(ys)>0){
    w <- max(ys$x)-min(ys$x)
    y <- (yu+yl)/2

    df0 <- data.frame(w=w,y=y)
    if(is.null(df)){df = df0}else{df=rbind(df,df0)}
  }}

wmax <- max(df$w)
ywidest <- mean(subset(df, w %in% wmax)$y)
wmin <- min(subset(df, y < ywidest)$w)
ynarrowest <-  mean(subset(df, w %in% wmin)$y)
ggplot()+
  geom_polygon(data= shape, aes(y=x+0.5,x=y,group=paste(part,hole), fill=paste(part,hole)))+
  geom_line(data= df, aes(x=y,y=w))+
  scale_x_continuous(breaks = c(0:5)/5)+
  scale_y_continuous(breaks = c(0:5)/5)+
  coord_flip(xlim=c(0,1),ylim=c(0,1))




shape <- as.polygons(maple)
coorwidth <- shape |> geom() |> as.data.frame()
coorwidth <- max((max(coorwidth$x) - min(coorwidth$x)),(max(coorwidth$x) - min(coorwidth$x)))
shape <- simplifyGeom(shape, coorwidth/1000) |> geom() |> as.data.frame()
shape <- shape |> group_by(part) |> mutate(size = length(part))
maxpart <- max(shape$size)
shape <- subset(shape, hole <=0)
shape <- subset(shape, size %in% maxpart)
# shape <- shape |> geom() |> as.data.frame()

make_tree <- function(ht.max, ht.min, crwd, dbh, crshape, stshape){
  crown <- subset(shapes, shape %in% crshape) |> mutate(x=x*crwd, z=z*(ht.max-ht.min)+ht.min, obj='crown')
  base <- subset(shapes, shape %in% stshape) |> mutate(x=x*dbh/100*1.1, z=z*(ht.min), obj='stem')
  tree = rbind(crown, base)
  tree$ptord <- rownames(tree) |> as.numeric()
  return(tree)}


ggplot()+
  geom_polygon(data= shape, aes(y=x+0.5,x=y), fill='red')+
  # geom_point(data= shape, aes(y=x+0.5,x=y))+
  geom_line(data= df, aes(x=y,y=w))+
  scale_x_continuous(breaks = c(0:5)/5)+
  scale_y_continuous(breaks = c(0:5)/5)+
  coord_flip(xlim=c(0,1),ylim=c(0,1))


ht.max = 30
ht.min = 12
dbh = 15
crwd = 8
crshape = branch1

morph_tree <- function(ht.max, ht.min, crwd, dbh, crshape, stshape) {
  shape <- crshape
  shape$x <- (shape$x - (max(shape$x) + min(shape$x))/2) / (max(shape$x) - min(shape$x))
  shape$y <- 1-1*(shape$y - min(shape$y)) / (max(shape$y) - min(shape$y))
  shape$y <- shape$y*-1+1
  #first pass to get wide and narrow
  df=NULL
  for(i in 1:25){#i=10
    yu = (i+3)/25
    yl = (i-3)/25
    ys <- subset(shape, y>= yl & y<= yu)
    if(nrow(ys)>0){
      w <- max(ys$x)-min(ys$x)
      y <- (yu+yl)/2

      df0 <- data.frame(w=w,y=y )
      if(is.null(df)){df = df0}else{df=rbind(df,df0)}
    }}

  wmax <- max(df$w)
  ywidest <- mean(subset(df, w %in% wmax)$y)
  wmin <- min(subset(df, y < ywidest)$w)
  ynarrowest <-  mean(subset(df, w %in% wmin)$y)

 #second pass for breast height
  df=NULL
  for(i in 1:50){#i=10

    yu = (i+3)/50
    yl = (i)/50
    ys <- subset(shape, y>= yl & y<= yu)
    if(nrow(ys)>0){
      w <- max(ys$x)-min(ys$x)
      y <- yl

      df0 <- data.frame(w=w,y=y)
      if(is.null(df)){df = df0}else{df=rbind(df,df0)}
    }}

  cht <- ((ywidest-ynarrowest)*0.5+ynarrowest)
  bh <- cht*1.4/ht.min

  bhw <- min(subset(df, y <= bh, select=w))

  dbhincrease <- (dbh/100)/bhw

  crownincrease <- crwd/1

  bottomincrease <- ht.min/cht

  topincrease <- (ht.max-ht.min)/(1-cht)

  shapenew <- shape |> mutate(transition = (y-ynarrowest+0.0000001)/(cht-ynarrowest+0.0000001),
                              xn = case_when(y < ynarrowest ~ x*dbhincrease,
                                            y >  cht ~ x*crownincrease,
                                            TRUE ~ x*((1-transition)*dbhincrease + (transition)*crownincrease)),
                              yn = case_when(y <= cht ~ y*bottomincrease,
                                            y > cht ~ cht*bottomincrease+(y-cht)*topincrease) )

  return(shapenew)}

ht.max = 30
ht.min = 10
dbh = 100
crwd = 20
crshape = branch1
newtree <- morph_tree(ht.max=ht.max, ht.min=ht.min, crwd=crwd, dbh=dbh, crshape=crshape)
branch1 <- branch1 |> mutate(id=as.numeric(rownames(branch1)))

ggplot()+
  geom_polygon(data= newtree, aes(y=yn,x=xn), fill='red')+
  geom_point(data= newtree, aes(y=yn,x=xn))+
  geom_text(data= newtree, aes(y=yn,x=xn, label=id))+
  scale_x_continuous(breaks = c(-5:5)*5)+
  scale_y_continuous(breaks = c(-5:5)*5)+
  coord_fixed(xlim=c(-20,20),ylim=c(0,ht.max+0))
ggplot()+
  geom_polygon(data= branch1, aes(y=y,x=x), fill='red')+
  geom_point(data= branch1, aes(y=y,x=x))+
  geom_text(data= branch1, aes(y=y,x=x, label=id))+
  scale_x_continuous(breaks = c(-5:5)/10)+
  scale_y_continuous(breaks = c(-5:5)/10)+
  coord_fixed(xlim=c(-0.5,0.5),ylim=c(0,1))

