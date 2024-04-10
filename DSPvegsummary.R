library(vegnasis)
library(data.table)
veg0 <- vegnasis::nasis.veg
veg0 <- subset(veg0, vegplotname %in% c('C-RF-1','C-RF-2','C-RF-3'))

veg <- veg0 |> clean.veg() |> fill.type.df() |> fill.hts.df()
veg.c <- veg |> group_by(label, taxon) |> summarise(cover = sum(cover), ht.max = max(ht.max)) |> group_by(label) |> mutate(tcover = sum(cover)) |> ungroup() |> mutate(rcover = cover/tcover) 

veg.c <- veg.c |> arrange(label, -rcover)
initlabel <- veg.c$label[1]
veg.c <- veg.c |> mutate(ccover0=0,ccover=0)
for(i in 1:nrow(veg.c)){#i=1
  initlabel0 <- veg.c$label[i]
if(i==1 | !initlabel %in%  initlabel0){
  veg.c$ccover0[i] <- 0
  cov <- 0+veg.c$rcover[i]
  veg.c$ccover[i] <- cov
  initlabel <- initlabel0
}else{
  veg.c$ccover0[i] <- cov
  cov = cov+veg.c$rcover[i]
  veg.c$ccover[i] <- cov
  initlabel <- initlabel0
}}

veg.c <- subset(veg.c, ccover0 < 0.5)

veg.ass <- get.assoc(veg)

veg.str <- get.structure(veg) 



veg.str <- veg.str |> summarise(tree.min = min(tree), tree.rv = mean(tree), tree.max = max(tree),
                             shrub.min = min(shrub), shrub.rv = mean(shrub), shrub.max = max(shrub),
                             herb.min = min(herb), herb.rv = mean(herb), herb.max = max(herb),
                             moss.min = min(tree), moss.rv = mean(moss), moss.max = max(moss),
                             ht.min = min(ht.max), ht.rv = mean(ht.max), ht.max = max(ht.max))

veg.str <- t(veg.str)
BA.to.SI(10)
veg <- veg |> mutate(BA = case_when(
  is.na(BA) & label %in% 'C-RF-1' & !is.na(dbh.max) & taxon %in% c('Acer nigrum', 'Carya ovata','Ulmus americana') ~ BA.to.SI(10),
  is.na(BA) & label %in% 'C-RF-1'  & !is.na(dbh.max) & taxon %in% 'Fagus grandifolia' ~ BA.to.SI(30),
  TRUE ~ BA))
veg.d <- veg |> subset(!is.na(BA)) |> mutate(dbh = ifelse(is.na(dbh.min), dbh.max,(dbh.min+dbh.max)/2)) |> group_by(label) |> 
  summarise(qdbh = weighted.mean(dbh, BA, na.rm=T), BA=sum(BA, na.rm=T))

veg.d <- veg.d |> summarise(qdbh.min = min(qdbh), qdbh.rv = mean(qdbh), qdbh.max = max(qdbh),
                                BA.min = min(BA), BA.rv = mean(BA), BA.max = max(BA))

veg.summ <-  summary.ESIS(veg, breaks = c(5), lowerQ = 0, upperQ = 1) #Mutliple relevé summary.
write.csv(veg.summ, 'DSPveg.summ.csv', row.names = F)