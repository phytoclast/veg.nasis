#FS legacy data ----
#install and load package
remotes::install_github("phytoclast/vegnasis", dependencies = FALSE)
library(vegnasis)


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

veg <- rbind(veg1,veg2)
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

