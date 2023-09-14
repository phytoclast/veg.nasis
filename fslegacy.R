#FS legacy data ----
#install and load package
# remotes::install_github("phytoclast/vegnasis", dependencies = FALSE)

library(vegnasis)
library(vegan)
library(proxy)


#provisional function to convert cover classes to absolute cover.
ECS.class.convert <-  function(class){
  cover = case_when(class %in% 'R' ~ 1,
                    class %in% 'C' ~ 5,
                    class %in% 'A' ~ 25)}

#Set active directory to be the same as where script file resides
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load FS legacy data
leg <- read.csv('data/PlantsDataCleaned.csv')

#Begin to assign column names that apply to all strata.
leg <- leg |> mutate(plot = vegplotid,
                     symbol = PLANTS_CODE,
                     taxon = fill.taxon.from.symbols(symbol))

#Subset data into strata, and provide sources for additional standardized column names.
leg.over <- leg |> subset(OVERSTORY %in% c('R','C','A')) |> mutate(cover = ECS.class.convert(OVERSTORY),stratum.min = 5, stratum.max = 30)
leg.under <- leg |> subset(UNDERSTORY %in% c('R','C','A')) |> mutate(cover = ECS.class.convert(UNDERSTORY),stratum.min = 0.5, stratum.max = 5)
leg.ground <- leg |> subset(GROUNDCOVER %in% c('R','C','A')) |> mutate(cover = ECS.class.convert(GROUNDCOVER),stratum.min = 0, stratum.max = 0.5)

#Combine rows of data from each stratum subset then fill in standardized columns.
veg1 <-  rbind(leg.over,leg.under,leg.ground) |> pre.fill.veg() |>
  fill.type.df() |> fill.hts.df() |> mutate(taxon = harmonize.taxa(taxon, fix = TRUE)) |> fill.nativity.df()



veg2 <- vegnasis::nasis.veg |> clean.veg() |>
  subset(grepl('MI165', plot))|>
  fill.type.df() |> fill.hts.df() |> mutate(taxon = harmonize.taxa(taxon, fix = TRUE)) |> fill.nativity.df()

veg <- bind_rows(veg1,veg2)
veg.m <- make.plot.matrix(veg, 'log')
d = vegan::vegdist(veg.m, method='bray')
t <- cluster::agnes(d, method = 'ward')|> as.hclust()
ape::plot.phylo(ape::as.phylo(t))
t <- optpart::flexbeta(d, beta = -0.2)|> as.hclust()
ape::plot.phylo(ape::as.phylo(t))


k = 3
groups <- cutree(t, k = k)
groups <- dendrogrouporder(t, groups)
a <- 'Wexford Legacy Data'
export.dendro(a,d,t,groups)


veg <- veg |> mutate(genus = link.taxonomy(taxon, taxrank=1), family = link.taxonomy(taxon, taxrank=2))
veg.matrix1 <- veg |> make.plot.matrix(tr='log', rc=T, nr=T, taxon='taxon')
veg.matrix2 <- veg |> make.plot.matrix(tr='log', rc=T, nr=T, taxon='genus')
veg.matrix3 <- veg |> make.plot.matrix(tr='log', rc=T, nr=T, taxon='family')

newnames <- subset(veg, grepl('.PL$', plot), select = plot) |> unique()
newnames$braytaxon <- 0
newnames$braygenus <- 0
newnames$brayfamily <- 0
newnames$simptaxon <- 0
newnames$simpgenus <- 0
newnames$simpfamily <- 0
for(i in 1:nrow(newnames)){
thisname <- substr(newnames$plot[i] , 1,14)

veg.1 <- veg.matrix1 |> subset(grepl(thisname, rownames(veg.matrix1)))
veg.2 <- veg.matrix2 |> subset(grepl(thisname, rownames(veg.matrix2)))
veg.3 <- veg.matrix3 |> subset(grepl(thisname, rownames(veg.matrix3)))

newnames[i,]$braytaxon <- vegdist(veg.1, method='bray')[1]
newnames[i,]$simptaxon <- as.dist(simil(veg.1, method='Simpson'))[1]
newnames[i,]$braygenus <- vegdist(veg.2, method='bray')[1]
newnames[i,]$simpgenus <- as.dist(simil(veg.2, method='Simpson'))[1]
newnames[i,]$brayfamily <- vegdist(veg.3, method='bray')[1]
newnames[i,]$simpfamily <- as.dist(simil(veg.3, method='Simpson'))[1]
}

write.csv(newnames, 'newnames.csv', row.names = F)




rankplots <- function(d){
dtab <- as.matrix(d) |> as.data.frame()
dtab$plot <- rownames(dtab)
dtab <- subset(dtab, grepl('P$', plot), select = grepl('PL$',colnames(dtab)))
pnames <-  colnames(dtab)
newtab = NULL
for (i in 1:length(colnames(dtab))){#i=1
  thisname <- pnames[i]
thistab <- subset(dtab, select = c(thisname))
thistab$rnk <- 1-(rank(thistab[,1])-1)/nrow(thistab)

diss <- subset(thistab, rownames(thistab)  %in% substr(thisname , 1,14))[1,1]
prank <- subset(thistab, rownames(thistab)  %in% substr(thisname , 1,14))[1,'rnk']
newtab0 <- data.frame(plot=thisname, diss=diss, prank = prank)
if(is.null(newtab)){newtab <- newtab0}else{newtab <- rbind(newtab, newtab0)}
}
return(newtab)
}

d = vegan::vegdist(veg.matrix1, method='bray')
bray.taxa <- rankplots(d)
write.csv(bray.taxa, 'bray.taxa.csv', row.names = FALSE)
d = vegan::vegdist(veg.matrix2, method='bray')
bray.taxa <- rankplots(d)
write.csv(bray.taxa, 'bray.gen.csv', row.names = FALSE)
d = vegan::vegdist(veg.matrix3, method='bray')
bray.taxa <- rankplots(d)
write.csv(bray.taxa, 'bray.fam.csv', row.names = FALSE)


d = as.dist(proxy::simil(veg.matrix1, method='Simpson'))
sim.taxa <- rankplots(d)
write.csv(sim.taxa, 'sim.taxa.csv', row.names = FALSE)
d = as.dist(proxy::simil(veg.matrix2, method='Simpson'))
sim.taxa <- rankplots(d)
write.csv(sim.taxa, 'sim.gen.csv', row.names = FALSE)
d = as.dist(proxy::simil(veg.matrix3, method='Simpson'))
sim.taxa <- rankplots(d)
write.csv(sim.taxa, 'sim.fam.csv', row.names = FALSE)



