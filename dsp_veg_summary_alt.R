library(vegnasis)
library(climatools)
library(soilDB)
library(aqp)
library(tidyr)

condRound10 <- function(x){
  x <- ifelse(x < 0.05, 0, ifelse(x < 10, round(x, 1), round(x,0)))
  x <- as.character(x)
  return(x)
}
condRound1 <- function(x){
  x <-  ifelse(x < 0.05, 0, ifelse(x < 1, round(x, 1), round(x,0)))
  x <- as.character(x)
  return(x)
}

#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# veg.raw0 <- soilDB::get_vegplot_species_from_NASIS_db(SS=F)
# saveRDS(veg.raw0,'jackpine/veg.raw0.RDS')

# landuse <- read.csv('jackpine/dsp_landuse.csv')
# landuse$landuse <- str_replace_all(landuse$landuse, '-','_')
# veg.raw <- readRDS('jackpine/veg.raw0.RDS')
siteass <- get_site_association_from_NASIS(SS=F)
sites <- get_site_data_from_NASIS_db(SS=F)
veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS=F)
vegplot <- soilDB::get_vegplot_from_NASIS_db(SS=F)
# landuse <- data.frame(plot=unique(vegplot$vegplotid), landuse = 'Wet Prairie')
landuse <- data.frame(siteobsiid=sites$siteobsiid, landuse = sites$commphasename)
landuse <- landuse |> left_join(data.frame(siteobsiid=vegplot$siteobsiid, plot = vegplot$vegplotid)) |> subset(!is.na(plot) & !is.na(landuse))
veg <- clean.veg(veg.raw)|> subset(!is.na(taxon)) 
veg <- veg |> inner_join(landuse)

veg <- veg |> mutate(type=NA) |> fill.type.df() |> fill.hts.df()
veg <- veg |> mutate(taxon = harmonize.taxa(veg$taxon, fix = TRUE, sensu = "usda"))

#Get vegetation Structure ----
veg.str <- veg |> get.structure(simple = TRUE)



veg.str <- veg.str |> inner_join(landuse)
veg.str.long <- tidyr::pivot_longer(veg.str, c(tree,shrub,herb,moss,ht.max))
veg.str.summary <- veg.str.long |> group_by(landuse, name) |> summarise(Low = round(quantile(value,0.05),1),
                                                                 RV = round(mean(value),1),
                                                                 High = round(quantile(value,0.95),1))

veg.str.wide <- tidyr::pivot_wider(veg.str.summary, names_from = landuse, values_from = c(Low,RV,High)) |> as.data.frame()

library(kableExtra)
library(knitr)

df2 <- veg.str.summary |> group_by(landuse) |> arrange()
df2$name <- factor(df2$name, levels = c('ht.max','tree','shrub','herb','moss'))
df2 <- df2[order(df2$landuse,df2$name),]
df2 |>
  knitr::kable(row.names = FALSE, digits = c(0,0,0,0)) |>
  remove_column(1) |>
  kableExtra::group_rows(index = table(df2$landuse)) |>
  # kableExtra::kable_paper("hover", full_width = F)
kable_classic(full_width = F, html_font = "Cambria")

  


# write.csv(df2, 'jackpine/maumeeveg.str.wide.csv', row.names = F, na='')


#Species_Composition
taxon.fill <- merge(data.frame(group = unique(veg$landuse)), data.frame(taxon = unique(veg$taxon), Low = 0, RV = 0, High = 0)) |> mutate(type = vegnasis::fill.type(taxon)) |> unique() 
taxon.fill <- taxon.fill[,c('group','taxon', 'type', 'Low', 'RV', 'High')]

veg.comp.summary <-  veg  |> summary.ESIS(group='landuse', breaks = c(5), normalize = F,  
                                          lowerQ = 0, upperQ = 1) |> ungroup()
veg.comp.summary <- veg.comp.summary |> mutate(Low = cover.Low, RV=cover.mean, High=cover.High)
overstory <- veg.comp.summary |> subset(Top > 5, select = c("group","taxon", "type","Low","RV","High"))
#add missing rows
o2 <- subset(taxon.fill, taxon %in% overstory$taxon)
o2 <- subset(o2, !paste(taxon,group) %in% paste(overstory$taxon,overstory$group) )
overstory <- overstory |> rbind(o2)


df2 <- overstory 
allplots <- df2 |> group_by(taxon, type) |> summarise(group = "All Landuses", Low = min(Low), RV = mean(RV), High = max(High)) |> arrange(-RV )
factorgroup <- unique(df2$group)
factortaxon <- allplots$taxon
df2 <- rbind(df2, allplots)
df2$taxon <- factor(df2$taxon, levels = factortaxon)
df2$group <- factor(df2$group, levels = c(factorgroup,"All Landuses"))
df2 <- df2 |> arrange(group, taxon)

df2 |> #mutate(Low = condRound1(Low), RV = condRound1(RV), High = condRound1(High)) |>
  knitr::kable(row.names = FALSE, digits = c(1,1,1,1,1)) %>%
  remove_column(1) |> column_spec(1,italic=T) |>
  kableExtra::group_rows(index = table(df2$group)) |>
  # kableExtra::kable_paper("hover", full_width = F)
  kable_classic(full_width = F, html_font = "Cambria")




# write.csv(overstory, 'jackpine/overstory.csv', row.names = F, na='')

# veg <- veg |> mutate(family = vegnasis::link.taxonomy(taxon, taxrank=2))

understory <- veg.comp.summary |> subset(Top <= 5, select = c("group","taxon", "type", "Low","RV","High")) 
u2 <- subset(taxon.fill, taxon %in% understory$taxon)
u2 <- subset(u2, !paste(taxon,group) %in% paste(understory$taxon,understory$group))
understory <- understory |> rbind(u2)



df2 <- understory 
allplots <- df2 |> group_by(taxon, type) |> summarise(group = "All Landuses", Low = min(Low), RV = mean(RV), High = max(High)) |> arrange(-RV )
keeptaxa <- subset(allplots, High >= 10)$taxon
factorgroup <- unique(df2$group)
factortaxon <- allplots$taxon
df2 <- rbind(df2, allplots)
df2$taxon <- factor(df2$taxon, levels = factortaxon)
df2$group <- factor(df2$group, levels = c(factorgroup,"All Landuses"))
df2 <- df2 |> arrange(group, taxon) |> subset(taxon %in% keeptaxa)

options(knitr.kable.NA = '-')
df2 |> #mutate(Low = condRound1(Low), RV = condRound1(RV), High = condRound1(High)) |>
  knitr::kable(row.names = FALSE, digits = c(1,1,1,1,1)) %>% 
  remove_column(1) |> column_spec(1,italic=T) |>
  kableExtra::group_rows(index = table(df2$group)) |>
  kable_styling(bootstrap_options = 'bordered') |>
  column_spec(c(1:5), bold = TRUE, border_left = TRUE, border_right = TRUE, color = "black", background = "white")|> 
  # kableExtra::kable_paper("hover", full_width = F)
  kable_classic(full_width = F, html_font = "Cambria")



# write.csv(understory, 'jackpine/maumeeunderstory.csv', row.names = F, na='')
