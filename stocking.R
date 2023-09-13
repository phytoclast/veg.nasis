library(vegnasis)
# Dey, D.C., Kabrick, J.M. and Schweitzer, C.J., 2017. Silviculture to restore oak savannas and woodlands. Journal of Forestry, 115(3), pp.202-211.
ba = c(38.5,46.5,52.5,56.5,60,62.5,65,67.5,70,72.5)
dbh = c(3,4,5,6,7,8,9,10,12,14)
stocking = c(100,100,100,100,100,100,100,100,100,100)

tabstock = data.frame(ba=ba,dbh=dbh,stocking=stocking)
tabstock = tabstock |> mutate(ba = BA.to.SI(ba), dbh = dbh.metric(dbh))


ba = c(15,65.5,75,80,70,
       30,65,48.5,40,24,
       23,24.5,40)
dbh = c(16,14,16,25,20.5,
        16,28,30,24,24,
        14,30,26)
cover = c(20,90,100,100,90,
          40,80,60,50,30,
          40,30,50)

tabcover = data.frame(ba=ba,dbh=dbh,cover=cover)
tabcover = tabcover |> mutate(ba = BA.to.SI(ba), dbh = dbh.metric(dbh))


mod <- lm(ba ~ dbh+I(dbh^0.5), data=tabstock)
summary(mod)
get.stocking = function(dbh, BA){
 stock = 100*BA/(-7.5236 + -0.5875*dbh + 7.5331*dbh^0.5)
 return(stock)
}

get.stocking(dbh.metric(14),BA.to.SI(40))

library(ggplot2)
ggplot(aes(ba,dbh), data=tabcover)+
  geom_point()



mod <- lm(cover ~ I(ba^1)+I(dbh^0.5), data=tabcover)
summary(mod)


mod <- lm(cover ~ ba+I(1/dbh^2), data=tabcover)
summary(mod)

get.cover = function(dbh, BA){
  cover = pmin(-4.9343 + 5.5347*BA + 15507.2062/dbh^2, 100)
  return(cover)
}

get.cover(dbh.metric(20),BA.to.SI(30))

