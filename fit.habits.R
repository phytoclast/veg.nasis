library(stringr)
library(dplyr)
library(ranger)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



f1 <- read.delim('data/plants/FinalForms.txt')
f2 <- read.delim('data/plants/List_Species2011.txt')
f3 <- read.delim('data/plants/BinomialGrowthHabits.txt')
bm.geo <- read.csv('data/plants/bm.geo.csv')
syns <- read.csv('data/plants/m.ac.csv')


sort(unique(bm.geo$STATECODE))


FC.1 <- c("US.HI")
FC.2 <- c("MX.AG","MX.BN", "MX.BS", "MX.CA","MX.CH","MX.DF", "MX.DU", "MX.GJ","MX.NL","MX.SO", "MX.TL","MX.ZA","US.AZ","US.NM")
FC.3 <- c("MX.CL", "MX.CP","MX.GR","MX.HI","MX.JA", "MX.MC","MX.MR","MX.MX", "MX.NA","MX.OA","MX.PU", "MX.QE", "MX.SI","MX.SL", "MX.TM", "MX.VE")
FC.4 <- c("MX.CM", "MX.QR", "MX.TB", "MX.YU", "US.PR","US.UM","US.VI")
FC.5 <- c("CA.AB","CA.BC","CA.MB","CA.SK","US.CA", "US.CO","US.ID","US.MT","US.ND","US.NE","US.OR","US.NV","US.SD", "US.UT","US.WA","US.WY")
FC.6 <- c("CA.LB","CA.NB", "CA.NF", "CA.PE","CA.QU","CA.NS","CA.NT","CA.NU","CA.YT","GL","PM", "US.AK")
FC.7 <- c("US.AL", "US.AR", "US.FL", "US.GA","US.KS", "US.LA", "US.MS", "US.OK","US.TX","US.SC")
FC.8 <- c("CA.ON", "US.CT", "US.DE", "US.DC","US.IA", "US.IL", "US.IN", "US.KY","US.MA", "US.MD", "US.ME", "US.MI", "US.MN", "US.MO","US.NC","US.NH", "US.NJ",   "US.NY", "US.OH",   "US.PA", "US.RI", "US.TN", "US.WI", "US.VT","US.VA", "US.WV")

FC <- cbind(FC= "FC.1", STATECODE = FC.1 ) %>% 
  rbind(cbind(FC= "FC.2", FC.2 ))%>%
  rbind(cbind(FC= "FC.3", FC.3 ))%>%
  rbind(cbind(FC= "FC.4", FC.4 ))%>%
  rbind(cbind(FC= "FC.5", FC.5 ))%>%
  rbind(cbind(FC= "FC.6", FC.6 ))%>%
  rbind(cbind(FC= "FC.7", FC.7 ))%>%
  rbind(cbind(FC= "FC.8", FC.8 ))%>% as.data.frame()

bm.geo.FC <- bm.geo %>% left_join(FC) %>% subset(Status %in% "N")
bm.geo.sums <- bm.geo.FC %>% group_by(syn, FC) %>% summarise(ct = length(FC))
bm.geo.max <- bm.geo.sums %>% group_by(FC) %>% summarise(ctmax = max(ct))
bm.geo.sums <- bm.geo.sums %>% left_join(bm.geo.max) %>% mutate(pct = ct/ctmax*100)


FC.1.taxa <- subset(bm.geo, STATECODE %in% FC.1 & Status %in% "N", select = syn)[,1]

f1 <- f1 %>% mutate(FC1 = ifelse(Scientific.Name %in% FC.1.taxa, 1,0))


f1$FinalHabits <- as.factor(f1$FinalHabits)
f1$Genus <- as.factor(f1$Genus)
f1$Family <- as.factor(f1$Family)
f1$Order <- as.factor(f1$Order)
f1$Superorder <- as.factor(f1$Superorder)
f1$Subclass <- as.factor(f1$Class)
f1$Class <- as.factor(f1$Class)
f1$Superclass <- as.factor(f1$Superclass)
f1$Subdivision <- as.factor(f1$Subdivision)
unique(colnames(f1))
rf <-  ranger(FinalHabits ~ Genus+Family+Order, data= f1)
f1$model <-  predictions(predict(rf, f1))