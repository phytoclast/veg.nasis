library(stringr)
library(dplyr)
library(ranger)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



f1 <- read.delim('data/plants/FinalForms.txt')
f2 <- read.delim('data/plants/List_Species2011.txt')
f3 <- read.delim('data/plants/BinomialGrowthHabits.txt')
bm.geo <- read.csv('data/plants/bm.geo.csv')
syns <- read.csv('data/plants/m.ac.csv')

#clean hydric ----
#Downloaded 2020 wetland plant list from Excel spreadsheet from Army Corps website, pasted into text file.
hydric <- read.delim('data/plants/hydric.taxa.txt')
hydric <- sapply(hydric, str_trim) |> as.data.frame()
colnames(hydric) <- colnames(hydric) |> stringr::str_replace('^X.','') |> stringr::str_replace('^\\.','') |> stringr::str_replace('\\.\\.$','')

replacehydric <- function(x){
  x = case_when(
    x %in% 'UPL' ~ 0,
    x %in% 'FACU' ~ 0.25,
    x %in% 'FAC' ~ 0.5,
    x %in% 'FACW' ~ 0.75,
    x %in% 'OBL' ~ 1,
    TRUE ~ NA_real_)}

hydric2 <- hydric |> mutate(across(2:11, replacehydric))
#Downloaded full PLANTS database as text file.
PLANTS <- read.csv('data/plants/PLANTSdownloadData.txt')
PLANTS.legit <- PLANTS |> subset(!grepl('auct.',Genera.Binomial.Author) &
                                   !grepl('illeg.',Genera.Binomial.Author) &
                                   !grepl(' non',Genera.Binomial.Author) &
                                   !grepl('auct',Trinomial.Author) &
                                   !grepl('illeg.',Trinomial.Author) &
                                   !grepl(' non',Trinomial.Author),
                                 select=c(Accepted.Symbol, Symbol, Scientific.Name))
PLANTS.legit1 <- PLANTS.legit |> mutate(binomial = paste(    str_split_fixed(Scientific.Name, '[[:blank:]]',3)[,1],
                                                            str_split_fixed(Scientific.Name, '[[:blank:]]',3)[,2]))

PLANTS.binomials <- PLANTS.legit1
PLANTS.binomials1 <- PLANTS.binomials
colnames(PLANTS.binomials1) <-  paste0(colnames(PLANTS.binomials1),".1")
PLANTS.binomials <- PLANTS.binomials |> left_join(PLANTS.binomials1, by=c('Accepted.Symbol'='Symbol.1'), multiple = "all")
hydric.syns <-  hydric2 |> left_join(PLANTS.binomials, by=c('Scientific.Name'='Scientific.Name'), multiple = "all")

missing <- PLANTS.binomials |> subset(!binomial.1 %in% hydric.syns$binomial.1)

library(vegnasis)
missing <- missing |> mutate(GH = get.habit.code(missing$Scientific.Name)) |> subset(!grepl('^N',GH) & str_count(binomial.1, '. .') >= 1)
#add missing as presumed UPL taxa
hydric3 <-  hydric2 |> bind_rows(data.frame(Scientific.Name = unique(missing$binomial.1), other=0))

write.csv(hydric3, 'data/plants/hydric.csv', row.names = F)
#write hydric function ----
# data('nasis.veg')
# x <- nasis.veg |> clean.veg()
# x <- x |> mutate(type=fill.type(taxon, type)) |> fill.hts.df()

# get.wetness <- function(x, region = 'NCNE'){
#   x <- x |> group_by(plot, taxon) |> summarise(cover=cover.agg(cover))
#
#   #Lookup default plant hydric indicator status ----
#   region = 'NCNE'
#   hydric <- hydric |> mutate(status = case_when(
#     region == 'AGCP' ~ AGCP,
#     region == 'AW' ~ AW,
#     region == 'CB' ~ CB,
#     region == 'EMP' ~ EMP,
#     region == 'HI' ~ HI,
#     region == 'MW' ~ MW,
#     region == 'NCNE' ~ NCNE,
#     region == 'WMVC' ~ WMVC,
#     region == 'AK' ~ AK,
#     TRUE ~ other))
#
#   hydric <- hydric |> mutate(
#     status = ifelse(is.na(status), rowMeans(select(hydric,
#                                                    c('AGCP','AW','CB','EMP','GP','HI','MW','NCNE','WMVC','AK','other')), na.rm = TRUE),status))
#
#
#   x$status0 = NA_real_
#   #first try straight join ----
#   x <- x |> left_join(hydric[,c('Scientific.Name','status')], by = c('taxon'='Scientific.Name'), multiple = 'first')
#   x <- x |> mutate(status0 = ifelse(is.na(status0)| status0 %in% NA_real_, status, status0))
#   x <- x[,1:4]
#   #then try synonym join ----
#   x <- x |> left_join(syns[,c('acc','syn')], by=c('taxon'='syn'), multiple = 'first') |> left_join(hydric[,c('Scientific.Name','status')], by = c('taxon'='Scientific.Name'), multiple = 'first')
#   x <- x |> mutate(status0 = ifelse(is.na(status0)| status0 %in% NA_real_, status, status0))
#   x <- x[,1:4]
#   #finally try genus only ----
#   #not implemented
#
#   x <- x|> subset(!is.na(status0) & !is.na(cover)) |> group_by(plot) |> summarise(wetness = sum(cover*status0+0.0005)/sum(cover+0.001))
#   return(x)
# }
#





#clean gho ----


gho <- read.delim('data/plants/GrowthHabitOptions.txt')
gho <-  gho %>%
  mutate(First = substr(Revised.Symbol, 1,1), Second = substr(Revised.Symbol, 2,2), Last = substr(Revised.Symbol, 3,5))


fsyn <- f1 |> left_join(syns[,c('acc','syn')], by= c('Scientific.Name'='syn')) |> mutate(Scientific.Name = ifelse(is.na(acc),Scientific.Name, acc))

f1 <-  fsyn %>% left_join(gho[,c("Revised.Symbol","First","Second","Last")], by=c("FinalHabits"="Revised.Symbol"))
ferns <- subset(f1, Subdivision %in% c('Lycopodiophytina','Polypodiophytina'), select = 'Genus') |> unique()

f1$First <- as.factor(f1$First)
f1$Second <- as.factor(f1$Second)
#f1$Third <- as.factor(f1$Third)
f1$Last <- as.factor(f1$Last)
f1$Genus <- as.factor(f1$Genus)
f1$Family <- as.factor(f1$Family)
f1$Order <- as.factor(f1$Order)
f1$Superorder <- as.factor(f1$Superorder)
f1$Subclass <- as.factor(f1$Class)
f1$Class <- as.factor(f1$Class)
f1$Superclass <- as.factor(f1$Superclass)
f1$Subdivision <- as.factor(f1$Subdivision)

gho <- gho |> mutate(ht.max = case_when(
  First %in%  "T" & Second %in% c('2','') ~ 24,
  First %in%  "T" & Second %in% c('1') ~ 12,
  First %in%  c("L","E")  ~ 12,
  First %in%  "S" & Second %in% c('2','') ~ 3,
  First %in%  "S" & Second %in% c('1') ~ 0.3,
  First %in%  "H" ~ 0.6,
  First %in%  "N" ~ 0,
    TRUE ~ 0))

FC.1 <- c("US.HI")
FC.2 <- c("MX.AG","MX.BN", "MX.BS", "MX.CA","MX.CH","MX.DF", "MX.DU", "MX.GJ","MX.NL","MX.SO", "MX.TL","MX.ZA","US.AZ","US.NM")
FC.3 <- c("MX.CL", "MX.CP","MX.GR","MX.HI","MX.JA", "MX.MC","MX.MR","MX.MX", "MX.NA","MX.OA","MX.PU", "MX.QE", "MX.SI","MX.SL", "MX.TM", "MX.VE")
FC.4 <- c("MX.CM", "MX.QR", "MX.TB", "MX.YU", "US.PR","US.UM","US.VI")
FC.5 <- c("CA.AB","CA.BC","CA.MB","CA.SK","US.CA", "US.CO","US.ID","US.MT","US.ND","US.NE","US.OR","US.NV","US.SD", "US.UT","US.WA","US.WY")
FC.6 <- c("CA.LB","CA.NB", "CA.NF", "CA.PE","CA.QU","CA.NS","CA.NT","CA.NU","CA.YT","GL","PM", "US.AK")
FC.7 <- c("US.AL", "US.AR", "US.FL", "US.GA","US.KS", "US.LA", "US.MS", "US.OK","US.TX","US.SC")
FC.8 <- c("CA.ON", "US.CT", "US.DE", "US.DC","US.IA", "US.IL", "US.IN", "US.KY","US.MA", "US.MD", "US.ME", "US.MI", "US.MN", "US.MO","US.NC","US.NH", "US.NJ",   "US.NY", "US.OH",   "US.PA", "US.RI", "US.TN", "US.WI", "US.VT","US.VA", "US.WV")
#native ranges ----
FC <- cbind(FC= "FC.1", STATECODE = FC.1 ) %>%
  rbind(cbind(FC= "FC.2", FC.2 ))%>%
  rbind(cbind(FC= "FC.3", FC.3 ))%>%
  rbind(cbind(FC= "FC.4", FC.4 ))%>%
  rbind(cbind(FC= "FC.5", FC.5 ))%>%
  rbind(cbind(FC= "FC.6", FC.6 ))%>%
  rbind(cbind(FC= "FC.7", FC.7 ))%>%
  rbind(cbind(FC= "FC.8", FC.8 ))%>% as.data.frame()

bm.geo.FC <- bm.geo %>% left_join(FC) %>% subset(Status %in% "N")
bm.geo.sums <- bm.geo.FC %>% group_by(ac.binomial, FC) %>% summarise(ct = length(FC))
bm.geo.max <- bm.geo.sums %>% group_by(ac.binomial) %>% summarise(ctsum = sum(ct))
bm.geo.sums <- bm.geo.sums %>% left_join(bm.geo.max) %>% mutate(pct = ct/ctsum*100)

bm.geo.col <-  as.data.frame(cbind(ac.binomial=unique(bm.geo.sums$ac.binomial))) %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.1=pct), FC %in% 'FC.1', select=c(ac.binomial, FC.1)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.2=pct), FC %in% 'FC.2', select=c(ac.binomial, FC.2)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.3=pct), FC %in% 'FC.3', select=c(ac.binomial, FC.3)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.4=pct), FC %in% 'FC.4', select=c(ac.binomial, FC.4)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.5=pct), FC %in% 'FC.5', select=c(ac.binomial, FC.5)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.6=pct), FC %in% 'FC.6', select=c(ac.binomial, FC.6)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.7=pct), FC %in% 'FC.7', select=c(ac.binomial, FC.7)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(FC.8=pct), FC %in% 'FC.8', select=c(ac.binomial, FC.8)))

f1.geo <- f1 %>% left_join(bm.geo.col, by=c('Binomial' = 'ac.binomial')) %>%
  mutate(FC.1 = ifelse(is.na(FC.1),0,FC.1),
         FC.2 = ifelse(is.na(FC.2),0,FC.2),
         FC.3 = ifelse(is.na(FC.3),0,FC.3),
         FC.4 = ifelse(is.na(FC.4),0,FC.4),
         FC.5 = ifelse(is.na(FC.5),0,FC.5),
         FC.6 = ifelse(is.na(FC.6),0,FC.6),
         FC.7 = ifelse(is.na(FC.7),0,FC.7),
         FC.8 = ifelse(is.na(FC.8),0,FC.8))
#exotic ranges ----

FC <- cbind(FC= "FC.1", STATECODE = FC.1 ) %>%
  rbind(cbind(FC= "xFC.2", FC.2 ))%>%
  rbind(cbind(FC= "xFC.3", FC.3 ))%>%
  rbind(cbind(FC= "xFC.4", FC.4 ))%>%
  rbind(cbind(FC= "xFC.5", FC.5 ))%>%
  rbind(cbind(FC= "xFC.6", FC.6 ))%>%
  rbind(cbind(FC= "xFC.7", FC.7 ))%>%
  rbind(cbind(FC= "xFC.8", FC.8 ))%>% as.data.frame()

bm.geo.FC <- bm.geo %>% left_join(FC)
bm.geo.sums <- bm.geo.FC %>% group_by(ac.binomial, FC) %>% summarise(ct = length(FC))
bm.geo.max <- bm.geo.sums %>% group_by(ac.binomial) %>% summarise(ctsum = sum(ct))
bm.geo.sums <- bm.geo.sums %>% left_join(bm.geo.max) %>% mutate(pct = ct/ctsum*100)

bm.geo.col <-  as.data.frame(cbind(ac.binomial=unique(bm.geo.sums$ac.binomial))) %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.1=pct), FC %in% 'xFC.1', select=c(ac.binomial, xFC.1)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.2=pct), FC %in% 'xFC.2', select=c(ac.binomial, xFC.2)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.3=pct), FC %in% 'xFC.3', select=c(ac.binomial, xFC.3)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.4=pct), FC %in% 'xFC.4', select=c(ac.binomial, xFC.4)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.5=pct), FC %in% 'xFC.5', select=c(ac.binomial, xFC.5)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.6=pct), FC %in% 'xFC.6', select=c(ac.binomial, xFC.6)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.7=pct), FC %in% 'xFC.7', select=c(ac.binomial, xFC.7)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(xFC.8=pct), FC %in% 'xFC.8', select=c(ac.binomial, xFC.8)))


f1.geo <- f1.geo %>% left_join(bm.geo.col, by=c('Binomial' = 'ac.binomial')) %>%
  mutate(xFC.1 = ifelse(is.na(xFC.1),0,xFC.1),
         xFC.2 = ifelse(is.na(xFC.2),0,xFC.2),
         xFC.3 = ifelse(is.na(xFC.3),0,xFC.3),
         xFC.4 = ifelse(is.na(xFC.4),0,xFC.4),
         xFC.5 = ifelse(is.na(xFC.5),0,xFC.5),
         xFC.6 = ifelse(is.na(xFC.6),0,xFC.6),
         xFC.7 = ifelse(is.na(xFC.7),0,xFC.7),
         xFC.8 = ifelse(is.na(xFC.8),0,xFC.8))



rf <-  ranger(First ~ Genus+Family+Order+Superorder+Subclass+Class+Superclass+Subdivision+
                FC.1+FC.2+FC.3+FC.4+FC.5+FC.6+FC.7+FC.8+
                xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8, always.split.variables = 'Genus', data= subset(f1.geo, !First %in% c("",".")),
              respect.unordered.factors = T)

f1.geo$mod1 <-  predictions(predict(rf, f1.geo))

f1.geo <- f1.geo %>% mutate(mod1 = ifelse(First %in% "", as.character(mod1), as.character(First)) %>% as.factor())

rf <-  ranger(Second ~ mod1+Genus+Family+Order+Superorder+Subclass+Class+Superclass+Subdivision+
                FC.1+FC.2+FC.3+FC.4+FC.5+FC.6+FC.7+FC.8+
                xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8, always.split.variables = c('Genus','mod1'), data= subset(f1.geo, !Second %in% c("",".")),
              respect.unordered.factors = T)

f1.geo$mod2 <-  predictions(predict(rf, f1.geo))


rf <-  ranger(Last ~ mod1+Genus+Family+Order+Superorder+Subclass+Class+Superclass+Subdivision+
                FC.1+FC.2+FC.3+FC.4+FC.5+FC.6+FC.7+FC.8+
                xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8, always.split.variables = c('Genus','mod1'), data= subset(f1.geo, !Last %in% c("","B")),
              respect.unordered.factors = T)

f1.geo$mod3 <-  predictions(predict(rf, f1.geo))

f1.geo <- f1.geo %>% mutate(mod2 = ifelse(mod1 %in% c('E','L','N'), ".",as.character(mod2)), GH = paste0(mod1, mod2, mod3))

taxon.habits <- subset(f1.geo, select = c(Scientific.Name, FinalHabits,Last, mod1,mod2,mod3,GH))

colnames(taxon.habits) <- c('Scientific.Name', 'preGH', 'Last','Stem','Size','Leaf', 'GH')
ghtest <- subset(taxon.habits, Last != Leaf & !Last %in% "")
write.csv(ghtest, 'data/plants/ghtest.csv', row.names = F)
# ghtest2 <- read.csv('data/plants/ghtest2.csv') %>% filter(Leafnew != Leaf)
taxon.habits <- taxon.habits |> mutate(genus = str_split_fixed(Scientific.Name , '[[:blank:]]',3)[,1])

taxon.habits <- taxon.habits |> mutate(Stem = ifelse(genus %in% c('Cuscuta', 'Dendrophthora', 'Phoradendron'), "E", as.character(Stem)),
                                       Size = ifelse(genus %in% c('Cuscuta', 'Dendrophthora', 'Phoradendron'), ".", as.character(Size)),
                                       Leaf = ifelse(!Last %in% c("","B","G"), as.character(Last), as.character(Leaf)),
                                       Leaf = ifelse(genus %in% c('Cuscuta', 'Dendrophthora', 'Phoradendron'), "i", as.character(Leaf)),
                                       Stem = ifelse(genus %in% c('Tillandsia', 'Epidendrum'), "E", as.character(Stem)),
                                       Size = ifelse(genus %in% c('Tillandsia', 'Epidendrum'), ".", as.character(Size)),
                                       Leaf = ifelse(genus %in% c('Tillandsia', 'Epidendrum'), "F", as.character(Leaf)),
                                       Leaf = ifelse(genus %in% c('Andropogon'), "GW", as.character(Leaf)),
                                       Leaf = ifelse(genus %in% ferns$Genus & Stem %in% 'H' & Leaf %in% 'F', "FE", as.character(Leaf)),
                                       Size = ifelse(Leaf %in% 'FE' & Stem %in% 'H', '.', as.character(Size)),


                                       GH = paste0(Stem, Size, Leaf))
taxon.habits <- subset(taxon.habits, select = c(-Last)) |> unique()

#bring in max heights ----

ht.round <- function(ht){
  ifelse(ht >= 8,round(ht,0),
         ifelse(ht >= 3, floor(ht*2+0.499)/2,round(ht,1)
         ))}

plant.hts <- read.delim('data/plants/Plant_heights.txt')
dbplanthts <- read.delim('data/plants/databaseherbhts.txt')
colnames(dbplanthts) <- colnames(plant.hts)
dbplanthts$Ht_m <- ht.round(dbplanthts$Ht_m)
dbplanthts <- subset(dbplanthts, !Scientific.Name %in% plant.hts$Scientific.Name & grepl(' ',Scientific.Name))
plant.hts <-rbind(plant.hts, dbplanthts)


x <- taxon.habits  |> mutate(ht.max = NA_real_)
x <- x |> left_join(plant.hts[,c('Scientific.Name','Ht_m')], by = c('Scientific.Name'='Scientific.Name'), multiple = 'first')
x <- x |> mutate(ht.max = Ht_m)
x <- x[,1:8]
x <- x |> left_join(syns[,c('acc','syn')], by=c('Scientific.Name'='syn'), multiple = 'first') |> left_join(plant.hts[,c('Scientific.Name','Ht_m')], by = c('acc'='Scientific.Name'), multiple = 'first')
x <- x |> mutate(ht.max = ifelse(is.na(ht.max), Ht_m, ht.max))
x <- x[,1:8]


x <- x  |> mutate(ht.max = case_when(
  # is.na(ht.max) & Stem %in%  "T" & Size %in% c('2','') ~ 24,
  # is.na(ht.max) & Stem %in%  "T" & Size %in% c('1') ~ 12,
  # is.na(ht.max) & Stem %in%  c("L","E")  ~ 12,
  # is.na(ht.max) & Stem %in%  "S" & Size %in% c('2','') ~ 3,
  # is.na(ht.max) & Stem %in%  "S" & Size %in% c('1') ~ 0.3,
  # is.na(ht.max) & Stem %in%  "H" & Leaf %in% "FV" ~ 3,
  # is.na(ht.max) & Stem %in%  "H" & Leaf %in% "A" ~ 0,
  # is.na(ht.max) & Stem %in%  "H" ~ 0.6,
  # is.na(ht.max) & Stem %in%  "N" ~ 0,
  ht.max <= 15 & Stem %in%  "T" & Size %in% c('2','') ~ 20,
  ht.max <= 5 & Stem %in%  "T" & Size %in% c('1') ~ 6,
  ht.max <= 0.5 & Stem %in%  "S" & Size %in% c('2','') ~ 1,
  ht.max > 15 & Stem %in%  "T" & Size %in% c('1') ~ 15,
  ht.max > 5 & Stem %in%  "S" & Size %in% c('2','') ~ 5,
  ht.max > 0.5 & Stem %in%  "S" & Size %in% c('1') ~ 0.5,
  TRUE ~ ht.max))
x <- x <- x  |> mutate(ht.train = case_when(
  is.na(ht.max) & Stem %in%  "T" & Size %in% c('2','') ~ 24,
  is.na(ht.max) & Stem %in%  "T" & Size %in% c('1') ~ 12,
  is.na(ht.max) & Stem %in%  c("L","E")  ~ 12,
  is.na(ht.max) & Stem %in%  "S" & Size %in% c('2','') ~ 3,
  is.na(ht.max) & Stem %in%  "S" & Size %in% c('1') ~ 0.3,
  is.na(ht.max) & Stem %in%  "H" & Leaf %in% "FV" ~ 3,
  is.na(ht.max) & Stem %in%  "H" & Leaf %in% "A" ~ 0,
  is.na(ht.max) & Stem %in%  "H" ~ 0.6,
  is.na(ht.max) & Stem %in%  "N" ~ 0,
  TRUE ~ ht.max),
  wt = ifelse(is.na(ht.max), 1,10000))

f1.geo2 <- f1.geo[,c(1:31)] |> left_join(x,multiple = "all") |> mutate(StemSize=paste0(Stem,Size))

rf <-  ranger(ht.train ~ StemSize+Stem+GH+Genus+Family+Order+Superorder+Subclass+Class+Superclass+Subdivision#+
                  # FC.1+FC.2+FC.3+FC.4+FC.5+FC.6+FC.7+FC.8+
                  # xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8
                , always.split.variables = c('Genus','Family','StemSize','GH','Stem'),mtry = NULL,
                data= subset(f1.geo2, !is.na(ht.train)), num.trees = 1000, sample.fraction = 0.1,
                respect.unordered.factors = T, case.weights = subset(f1.geo2, !is.na(ht.train), select = wt))

f1.geo2 <- f1.geo2 |> mutate(ht.max.pred = predictions(predict(rf, f1.geo2)))

x <- f1.geo2  |> mutate(ht.max = case_when(
  is.na(ht.max)  ~ ht.max.pred,
  is.na(ht.max) & Stem %in%  "T" & Size %in% c('2','') ~ 24,
  is.na(ht.max) & Stem %in%  "T" & Size %in% c('1') ~ 12,
  is.na(ht.max) & Stem %in%  c("L","E")  ~ 12,
  is.na(ht.max) & Stem %in%  "S" & Size %in% c('2','') ~ 3,
  is.na(ht.max) & Stem %in%  "S" & Size %in% c('1') ~ 0.3,
  is.na(ht.max) & Stem %in%  "H" & Leaf %in% "FV" ~ 3,
  is.na(ht.max) & Stem %in%  "H" & Leaf %in% "A" ~ 0,
  is.na(ht.max) & Stem %in%  "H" ~ 0.6,
  is.na(ht.max) & Stem %in%  "N" ~ 0,
  ht.max <= 15 & Stem %in%  "T" & Size %in% c('2','') ~ 20,
  ht.max <= 5 & Stem %in%  "T" & Size %in% c('1') ~ 6,
  ht.max <= 0.5 & Stem %in%  "S" & Size %in% c('2','') ~ 1,
  ht.max > 15 & Stem %in%  "T" & Size %in% c('1') ~ 15,
  ht.max > 5 & Stem %in%  "S" & Size %in% c('2','') ~ 5,
  ht.max > 0.5 & Stem %in%  "S" & Size %in% c('1') ~ 0.5,
  TRUE ~ ht.max),
  ht.max = ht.round(ht.max)) |> subset(select=c(Scientific.Name, ht.max)) |> unique()


# linmod <-  lm(ht.train ~ StemSize+Stem+GH+Genus+Family+Order+Superorder+Subclass+Class+Superclass+Subdivision+
#                 FC.1+FC.2+FC.3+FC.4+FC.5+FC.6+FC.7+FC.8+
#                 xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8
#               ,
#               data= subset(f1.geo2, !is.na(ht.train)), weights = subset(f1.geo2, !is.na(ht.train), select = wt)[,1])
#
# f1.geo2$ht.max.predlm <-  predict(linmod, f1.geo2, na.exclude = T)
#
#
# summary(linmod)


taxon.habits <- taxon.habits |>  unique() |> left_join(x)
















write.csv(taxon.habits, 'data/plants/taxon.habits.csv', row.names = F)
write.csv(gho, 'data/plants/gho.csv', row.names = F)

#append habits ----
taxon.habits <-  read.csv('data/plants/taxon.habits.csv')

bm.geo.FC <- bm.geo %>% left_join(FC) |> subset(!ac.binomial %in% "" & !is.na(ac.binomial))
bm.geo.sums <- bm.geo.FC %>% group_by(ac.binomial, FC) %>% summarise(ct = length(FC))
bm.geo.max <- bm.geo.sums %>% group_by(FC) %>% summarise(fcsum = max(ct))
bm.geo.sums <- bm.geo.sums %>% left_join(bm.geo.max) %>% mutate(pct = ct/fcsum*100)



taxon.habits.geo <- taxon.habits %>% left_join(bm.geo.sums, by=c('Scientific.Name' = 'ac.binomial'), multiple = "all") |> mutate(pct = ifelse(is.na(pct), 0.5,pct))

taxon.habits.genus.stem <- taxon.habits.geo |> group_by(genus, Stem) |> summarise(pct= sum(pct))
taxon.habits.genus.stem.max <- taxon.habits.genus.stem |> group_by(genus) |> summarise(pctmax= max(pct))
taxon.habits.genus.stem <- taxon.habits.genus.stem |> left_join(taxon.habits.genus.stem.max) |> subset(pctmax == pct)




taxon.habits.genus.GH <- taxon.habits.geo |> inner_join(taxon.habits.genus.stem[,c('genus','Stem')]) |> group_by(genus, GH) |> summarise(ht.max= weighted.mean(ht.max, pct), pct= sum(pct))
taxon.habits.genus.GH.max <- taxon.habits.genus.GH |> group_by(genus) |> summarise(pctmax= max(pct))
taxon.habits.genus.GH <- taxon.habits.genus.GH |> left_join(taxon.habits.genus.GH.max) |> subset(pctmax == pct)
taxon.habits.genus.GH$ht.max <- ht.round(taxon.habits.genus.GH$ht.max)




genus.habits <- f2 |> subset(!AcGenus %in% "" & grepl('^N',Form), select = c(AcGenus, Form)) |> unique()
genus.habits <- genus.habits |> mutate(genus = AcGenus, GH = ifelse(Form %in% 'NB', "N.B", "N.L"), ht.max= 0)
genus.habits <- genus.habits |> subset(select=c(genus, GH, ht.max)) |> rbind(taxon.habits.genus.GH[,c("genus","GH","ht.max")]) |>
  rbind(cbind(genus=c('Nostoc','Chara','Lyngbya','Agardhiella', 'Alaria', 'Ascophyllum','Bryopsis',
                      'Ceramium', 'Chaetomorpha', 'Champia', 'Chondrus', 'Cladophora', 'Codium',
                      'Colpomenia', 'Desmarestia', 'Ectocarpus', 'Enteromorpha', 'Fucus', 'Gracilaria',
                      'Hildenbrandtia', 'Hypnea', 'Laminaria','Nitella', 'Phyllitis', 'Polysiphonia',
                      'Porphyra', 'Rhizoclonium', 'Spyridia', 'Ulva', 'Vaucheria','Postelsia','Phyllospadix', 'Nemalion','Callithamnion','Pylaiella'), GH='N.A', ht.max=0))

write.csv(genus.habits, 'data/plants/genus.habits.csv', row.names = F)





















taxa <- c('Osmunda cinnamomea', 'Osmunda sp', 'Osmundastrum cinnamomeum', 'Washingtonia divaricata', 'Washingtonia filifera', 'Washingtonia', 'Washingtonia rubra', '')

get.habit.code <- function(taxa){
x  <-  as.data.frame(cbind(taxa=taxa)) |> mutate(GH0 = '', genus = str_split_fixed(taxa , '[[:blank:]]',3)[,1])
#first try straight join ----
x <- x |> left_join(taxon.habits[,c('Scientific.Name','GH')], by = c('taxa'='Scientific.Name'), multiple = 'first')
x <- x |> mutate(GH0 = ifelse(is.na(GH0)| GH0 %in% "", GH, as.character(GH0)))
x <- x[,1:3]
#then try synonym join ----
x <- x |> left_join(syns[,c('acc','syn')], by=c('taxa'='syn'), multiple = 'first') |> left_join(taxon.habits[,c('Scientific.Name','GH')], by = c('acc'='Scientific.Name'), multiple = 'first')
x <- x |> mutate(GH0 = ifelse(is.na(GH0)| GH0 %in% "", GH, as.character(GH0)))
x <- x[,1:3]
#finally try genus only ----
x <- x |> left_join(genus.habits, by = c('genus'='genus'), multiple = 'first')
x <- x |> mutate(GH0 = ifelse(is.na(GH0)| GH0 %in% "", GH, as.character(GH0)))
x <- x[,1:2]
return(x$GH0) }
code=c('T.BE','T2BD', 'E.i','E.F')
type = c('name', 'stem', 'ESIS')
get.habit <- function(code,type='name'){
  x  <-  as.data.frame(cbind(code=code)) |>
    left_join(gho, by = c('code'='Revised.Symbol'), multiple = 'first')
  if(type %in% 'stem'){return(x$First)}else if(type %in% 'ESIS'){return(x$ESIS.Group)}else{return(x$Habitname)}
}

get.habit(type='stem', code=code)
get.habit.code('Asimina triloba') |> get.habit('ESIS')


#Get nativity ----
#
#

Northwest <- c("US.WA","US.OR","US.ID","US.MT","US.WY")
Southwest <- c("US.NV","US.UT","US.CO","US.CA","US.AZ","US.NM")
NorthCentral <- c("US.IA","US.MO","US.MN","US.ND","US.NE","US.SD")
Southcentral <- c("US.AR","US.LA","US.KS", "US.OK","US.TX")
Northeast <- c("US.CT", "US.DE", "US.DC", "US.IL", "US.IN", "US.KY","US.MA", "US.MD", "US.ME", "US.MI","US.NH", "US.NJ","US.NY", "US.OH","US.PA", "US.RI",  "US.WI", "US.VT","US.VA", "US.WV")
Southeast <- c("US.TN","US.NC","US.AL",  "US.FL", "US.GA",  "US.MS", "US.SC")
Alaska <- c("US.AK")
Hawaii <- c("US.HI")
Caribbean <- c( "US.PR","US.UM","US.VI")
CanadaWest <- c("CA.AB","CA.BC","CA.MB","CA.SK")
CanadaEast <- c("CA.LB","CA.NB", "CA.NF", "CA.PE","CA.QU","CA.NS", "CA.ON","PM")
Arctic <- c("CA.NT","CA.NU","CA.YT","GL")
Mexico <- c("MX.CM", "MX.QR", "MX.TB", "MX.YU","MX.AG","MX.BN", "MX.BS", "MX.CA","MX.CH","MX.DF", "MX.DU", "MX.GJ","MX.NL","MX.SO", "MX.TL","MX.ZA","MX.CL", "MX.CP","MX.GR","MX.HI","MX.JA", "MX.MC","MX.MR","MX.MX", "MX.NA","MX.OA","MX.PU", "MX.QE", "MX.SI","MX.SL", "MX.TM", "MX.VE")
FC <- cbind(FC= "Northwest", STATECODE = Northwest ) %>%
  rbind(cbind(FC= "Southwest", Southwest ))%>%
  rbind(cbind(FC= "NorthCentral", NorthCentral ))%>%
  rbind(cbind(FC= "Southcentral", Southcentral ))%>%
  rbind(cbind(FC= "Northeast", Northeast ))%>%
  rbind(cbind(FC= "Southeast", Southeast ))%>%
  rbind(cbind(FC= "Alaska", Alaska ))%>%
  rbind(cbind(FC= "Hawaii", Hawaii ))%>%
  rbind(cbind(FC= "Caribbean", Caribbean ))%>%
  rbind(cbind(FC= "CanadaWest", CanadaWest ))%>%
  rbind(cbind(FC= "CanadaEast", CanadaEast ))%>%
  rbind(cbind(FC= "Arctic", Arctic ))%>%
  rbind(cbind(FC= "Mexico", Mexico ))%>% as.data.frame()

bm.geo.FC <- bm.geo %>% left_join(FC) %>% subset(Status %in% "N")
bm.geo.sums <- bm.geo.FC %>% group_by(ac.binomial, FC) %>% summarise(ct = length(FC))
bm.geo.max <- bm.geo.sums %>% group_by(ac.binomial) %>% summarise(ctsum = sum(ct))
bm.geo.sums <- bm.geo.sums %>% left_join(bm.geo.max) %>% mutate(pct = ct/ctsum*100)

bm.geo.col <-  as.data.frame(cbind(ac.binomial=unique(bm.geo.sums$ac.binomial))) %>%
  left_join(subset(bm.geo.sums %>% mutate(Northwest=pct), FC %in% 'Northwest', select=c(ac.binomial, Northwest)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Southwest=pct), FC %in% 'Southwest', select=c(ac.binomial, Southwest)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(NorthCentral=pct), FC %in% 'NorthCentral', select=c(ac.binomial, NorthCentral)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Southcentral=pct), FC %in% 'Southcentral', select=c(ac.binomial, Southcentral)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Northeast=pct), FC %in% 'Northeast', select=c(ac.binomial, Northeast)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Southeast=pct), FC %in% 'Southeast', select=c(ac.binomial, Southeast)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Alaska=pct), FC %in% 'Alaska', select=c(ac.binomial, Alaska)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Hawaii=pct), FC %in% 'Hawaii', select=c(ac.binomial, Hawaii)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Caribbean=pct), FC %in% 'Caribbean', select=c(ac.binomial, Caribbean)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(CanadaWest=pct), FC %in% 'CanadaWest', select=c(ac.binomial, CanadaWest)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(CanadaEast=pct), FC %in% 'CanadaEast', select=c(ac.binomial, CanadaEast)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Arctic=pct), FC %in% 'Arctic', select=c(ac.binomial, Arctic)))  %>%
  left_join(subset(bm.geo.sums %>% mutate(Mexico=pct), FC %in% 'Mexico', select=c(ac.binomial, Mexico)))
missing <- bm.geo |> subset(!ac.binomial %in% bm.geo.col$ac.binomial, select = ac.binomial) |> unique()
bm.geo.col <- bm.geo.col |> bind_rows(missing)

nativity <- bm.geo.col %>%
  mutate(Northwest = ifelse(is.na(Northwest),0,1),
         Southwest = ifelse(is.na(Southwest),0,1),
         NorthCentral = ifelse(is.na(NorthCentral),0,1),
         Southcentral = ifelse(is.na(Southcentral),0,1),
         Northeast = ifelse(is.na(Northeast),0,1),
         Southeast = ifelse(is.na(Southeast),0,1),
         Alaska = ifelse(is.na(Alaska),0,1),
         Hawaii = ifelse(is.na(Hawaii),0,1),
         Caribbean = ifelse(is.na(Caribbean),0,1),
         CanadaWest = ifelse(is.na(CanadaWest),0,1),
         CanadaEast = ifelse(is.na(CanadaEast),0,1),
         Arctic = ifelse(is.na(Arctic),0,1),
         Mexico = ifelse(is.na(Mexico),0,1))

write.csv(nativity, 'data/plants/nativity.csv', row.names = F)

#Fix Rubus----
taxon.habits <-  read.csv('data/plants/taxon.habits.csv')
fixrubus <-  read.csv('data/plants/fixrubus.csv')
fixrubus <-  fixrubus |> inner_join(taxon.habits, by=c('taxon'='Scientific.Name'))
fixrubus <-  fixrubus |> mutate(Stem = substr(correctedform, 1,1),
                                Size = substr(correctedform, 2,2),
                                Leaf = substr(correctedform, 3,5),
                                GH = correctedform,
                                correctedform=NULL)
fixrubus <- fixrubus |> mutate(ht.max = case_when(Stem %in% 'H' ~ 0.1,
                                                  Size %in% 1 ~ 0.3,
                                                  Size %in% 2 ~ 1))
taxon.habits1 <- taxon.habits |> subset(!Scientific.Name %in% fixrubus$taxon )
colnames(fixrubus)[colnames(fixrubus) %in% "taxon"] <- 'Scientific.Name'
taxon.habits1 <- rbind(taxon.habits1, fixrubus)
write.csv(taxon.habits1, 'data/plants/taxon.habits.csv', row.names = F)


