library(stringr)
library(dplyr)
library(ranger)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



f1 <- read.delim('data/plants/FinalForms.txt')
f2 <- read.delim('data/plants/List_Species2011.txt')
f3 <- read.delim('data/plants/BinomialGrowthHabits.txt')
bm.geo <- read.csv('data/plants/bm.geo.csv')
syns <- read.csv('data/plants/m.ac.csv')
gho <- read.delim('data/plants/GrowthHabitOptions.txt')
gho <-  gho %>% 
  mutate(First = substr(Revised.Symbol, 1,1), Second = substr(Revised.Symbol, 2,2), Last = substr(Revised.Symbol, 3,5))

f1 <-  f1 %>% left_join(gho[,c("Revised.Symbol","First","Second","Last")], by=c("FinalHabits"="Revised.Symbol"))


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
                xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8, always.split.variables = 'Genus', data= subset(f1.geo, !First %in% c("",".")))

f1.geo$mod1 <-  predictions(predict(rf, f1.geo))

f1.geo <- f1.geo %>% mutate(mod1 = ifelse(First %in% "", as.character(mod1), as.character(First)) %>% as.factor()) 

rf <-  ranger(Second ~ mod1+Genus+Family+Order+Superorder+Subclass+Class+Superclass+Subdivision+
                FC.1+FC.2+FC.3+FC.4+FC.5+FC.6+FC.7+FC.8+
                xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8, always.split.variables = c('Genus','mod1'), data= subset(f1.geo, !Second %in% c("",".")))

f1.geo$mod2 <-  predictions(predict(rf, f1.geo))


rf <-  ranger(Last ~ mod1+Genus+Family+Order+Superorder+Subclass+Class+Superclass+Subdivision+
                FC.1+FC.2+FC.3+FC.4+FC.5+FC.6+FC.7+FC.8+
                xFC.1+xFC.2+xFC.3+xFC.4+xFC.5+xFC.6+xFC.7+xFC.8, always.split.variables = c('Genus','mod1'), data= subset(f1.geo, !Last %in% c("","B")))

f1.geo$mod3 <-  predictions(predict(rf, f1.geo))
