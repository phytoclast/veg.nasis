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


FC.1 <- c("US.HI",)
FC.2 <- c("MX.AG","MX.BN", "MX.BS", "MX.CA","MX.CH","MX.DF", "MX.DU", "MX.GJ","MX.NL","MX.SO", "MX.TL","MX.ZA","US.AZ","US.NM")
FC.3 <- c("MX.CL", "MX.CP","MX.GR","MX.HI","MX.JA", "MX.MC","MX.MR","MX.MX", "MX.NA","MX.OA","MX.PU", "MX.QE", "MX.SI","MX.SL", "MX.TM", "MX.VE")
FC.4 <- c("MX.CM", "MX.QR", "MX.TB", "MX.YU", "US.PR","US.UM","US.VI")
FC.5 <- c("CA.AB","CA.BC","CA.MB","CA.SK","US.CA", "US.CO","US.ID","US.MT","US.ND","US.NE","US.OR","US.NV","US.SD", "US.UT","US.WA","US.WY")
FC.6 <- c("CA.LB","CA.NB", "CA.NF", "CA.PE","CA.QU","CA.NS","CA.NT","CA.NU","CA.YT","GL","PM", "US.AK")
FC.7 <- c("US.AL", "US.AR", "US.FL", "US.GA","US.KS", "US.LA", "US.MS", "US.OK","US.TX","US.SC")
FC.8 <- c("CA.ON", "US.CT", "US.DE", "US.DC","US.IA", "US.IL", "US.IN", "US.KY","US.MA", "US.MD", "US.ME", "US.MI", "US.MN", "US.MO","US.NC","US.NH", "US.NJ",   "US.NY", "US.OH",   "US.PA", "US.RI", "US.TN", "US.WI", "US.VT","US.VA", "US.WV")

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