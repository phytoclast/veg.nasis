library(ggplot2)
library(vegnasis)
library(terra)
library(sf)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

shapes <- vegnasis::shapes
plength = 100
#Establish a stand 50 by 20 m.

#rando stand map
#
cw <- 5
set.seed(0)
for(i in 1:50){
  n <- i*20
  
  for(k in 1:10){
    x <- runif(n, 0,100)
    y <- runif(n, 0,100)
    standmap <-  data.frame(x=x,y=y, cw=cw/2)
    standmap <-  standmap |> vect(geom=c("x", "y"), keepgeom=TRUE)
    stbuff <- buffer(standmap, width=cw/2) |> aggregate() |> crop(ext(0,100,0,100))
    cover <- expanse(stbuff)/10000*100
    covertab0 <- data.frame(n=n,cw=cw,cover=cover)
    if(i==1 & k==1){covertab <- covertab0}else{covertab <- rbind(covertab,covertab0)}
  }
}

covertab <- covertab |> group_by(n) |> summarise(cover = mean(cover))
covertab <- covertab |> mutate(ideal = (1-((1-(cw/2)^2*pi/10000)^n))*100, mxcover = pmin(100,n*(cw/2)^2*pi/10000*100), opt=(mxcover+ideal)/2)


#fixed

stand <- make_hex_stand(plength/100,1) 
set.seed(0)
for(i in 1:50){
  n <- i*20
  
  for(k in 1:5){
    standmap <- stand[sample(1:nrow(stand), size = n, replace = FALSE, prob = stand$wt),]
    standmap <-  data.frame(x=standmap$x,y=standmap$y)
    standmap <-  standmap |> vect(geom=c("x", "y"), keepgeom=TRUE)
    stbuff <- buffer(standmap, width=cw/2) |> aggregate() |> crop(ext(0,100,0,100))
    cover <- expanse(stbuff)/10000*100
    covertab0 <- data.frame(n=n,cw=cw,cover=cover)
    if(i==1 & k==1){covertab2 <- covertab0}else{covertab2 <- rbind(covertab2,covertab0)}
  }
}
covertab2 <- covertab2 |> group_by(n) |> summarise(cover = mean(cover))
stand <- make_hex_stand(plength/100,1) 
set.seed(0)
for(i in 1:50){
  n <- i*20
  
  for(k in 1:5){
    standmap <- stand[sample(1:nrow(stand), size = n, replace = FALSE, prob = stand$wt^2),]
    standmap <-  data.frame(x=standmap$x,y=standmap$y)
    standmap <-  standmap |> vect(geom=c("x", "y"), keepgeom=TRUE)
    stbuff <- buffer(standmap, width=cw/2) |> aggregate() |> crop(ext(0,100,0,100))
    cover <- expanse(stbuff)/10000*100
    covertab0 <- data.frame(n=n,cw=cw,cover=cover)
    if(i==1 & k==1){covertab3 <- covertab0}else{covertab3 <- rbind(covertab3,covertab0)}
  }
}
covertab3 <- covertab3 |> group_by(n) |> summarise(cover = mean(cover))

ggplot()+
  geom_line(data=covertab, aes(x=n,y=cover), color='red')+
  geom_line(data=covertab, aes(x=n,y=ideal), color='black')+
  geom_line(data=covertab, aes(x=n,y=opt), color='black')+
  geom_line(data=covertab, aes(x=n,y=mxcover), color='black')+
  geom_line(data=covertab2, aes(x=n,y=cover), color='blue')+
  geom_line(data=covertab3, aes(x=n,y=cover), color='green')#+
#scale_x_continuous(limits = c(0,500))


cw = 7
scf = cw/9 
stand <- make_hex_stand(plength/100/scf,1*scf) 
standmap <- stand[sample(1:nrow(stand), size = 150, replace = FALSE, prob = stand$wt^0.1),]#
standmap <-  data.frame(x=standmap$x,y=standmap$y)
standmap <-  standmap |> vect(geom=c("x", "y"), keepgeom=TRUE)
stbuff <- buffer(standmap, width=cw/2) |> aggregate() |> crop(ext(0,100,0,100))

ggplot()+
  geom_point(data=stand, aes(x=xp,y=yp, size=wt))+
  geom_sf(data=st_as_sf(stbuff), fill='red', alpha=0.5)

#-----------new method to grow stand based on shade
plength = 100
cw = 9
scf = cw/9 
stand <- make_hex_stand(plength/100/scf,1*scf)

# nudgx <- runif(nrow(stand), min=-0.25, max=0.25)
# nudgy <- runif(nrow(stand), min=-0.25, max=0.25)
# stand <- stand |> mutate(xp=xp+nudgx, yp=yp+nudgy)
seeds <- c(0,1,2,3,4)
cws <- c(6,7,8,9,10)
for(k in 1:5){
  cw = cws[k]
  set.seed(seeds[k])
  stand <- stand |> mutate(dbh = 0, cwd = 0, shd=0, d=NA, cd = 100)
  for(i in 1:1000){#i=1
    empty <- stand#[stand$dbh == 0,]
    empty <- empty[sample(1:nrow(empty),size=1, replace = FALSE, prob = empty$wt),]$stumpid
    stand <- stand |> mutate(dbh = ifelse(stumpid %in% empty, 30,dbh),
                             cwd = ifelse(stumpid %in% empty, cw,cwd))
    newtree <- subset(stand, stumpid %in% empty)
    newtree <- newtree |> mutate(cwd = pmax(cwd/1,pmin(cwd, cd*1)))
    stand <- stand |> mutate(d = ((newtree$yp-yp)^2+(newtree$xp-xp)^2)^0.5,
                             shd = ifelse(d < newtree$cwd/2, shd+0.1,shd),
                             cd = pmin(d-newtree$cwd/2, cd, na.rm = TRUE),
                             wt = ifelse(shd <= 0, 1,ifelse(shd <= 0.1, 1, 1)))
    totcov <- subset(stand, xp >= 10 & xp <= 90 & yp >= 10 & yp <= 90) |> 
      mutate(shd =  ifelse(shd >0,1,0), stem =  ifelse(dbh >0,1,0))
    
    cdf0 <- data.frame(crowns=sum(totcov$stem)/6400*10000*(cw/2)^2*pi/100, tcov = mean(totcov$shd)*100)
    if(i==1){cdf <- cdf0}else{cdf <- rbind(cdf, cdf0)}
  }
  if(k==1){cdfx <- cdf}else{cdfx <- rbind(cdfx, cdf)}
}
cdfoverlap <- cdfx

plength = 100
cw = 9
scf = cw/9 
stand <- make_hex_stand(plength/100/scf,1*scf)
seeds <- c(0,1,2,3,4)
cws <- c(6,7,8,9,10)
for(k in 1:5){
  cw = cws[k]
  set.seed(seeds[k])
  stand <- stand |> mutate(dbh = 0, cwd = 0, shd=0, d=NA, cd = 100)
  for(i in 1:1000){#i=1
    empty <- stand[stand$dbh == 0,]
    empty <- empty[sample(1:nrow(empty),size=1, replace = FALSE, prob = empty$wt),]$stumpid
    stand <- stand |> mutate(dbh = ifelse(stumpid %in% empty, 30,dbh),
                             cwd = ifelse(stumpid %in% empty, cw,cwd))
    newtree <- subset(stand, stumpid %in% empty)
    newtree <- newtree |> mutate(cwd = pmax(cwd/1,pmin(cwd, cd*1)))
    stand <- stand |> mutate(d = ((newtree$yp-yp)^2+(newtree$xp-xp)^2)^0.5,
                             shd = ifelse(d < newtree$cwd/2, shd+0.1,shd),
                             cd = pmin(d-newtree$cwd/2, cd, na.rm = TRUE),
                             wt = ifelse(shd <= 0, 1,ifelse(shd <= 0.1, 0.1, 0.01)))
    totcov <- subset(stand, xp >= 10 & xp <= 90 & yp >= 10 & yp <= 90) |> 
      mutate(shd =  ifelse(shd >0,1,0), stem =  ifelse(dbh >0,1,0))
    
    cdf0 <- data.frame(crowns=sum(totcov$stem)/6400*10000*(cw/2)^2*pi/100, tcov = mean(totcov$shd)*100)
    if(i==1){cdf <- cdf0}else{cdf <- rbind(cdf, cdf0)}
  }
  if(k==1){cdfx <- cdf}else{cdfx <- rbind(cdfx, cdf)}
}

cdfavoid <- cdfx
ggplot(cdfx, aes(x=crowns,y=tcov))+
  geom_point(size=0.1)



ccv <- function(x){
  #avoid overlap
  b1 = 0.002985  ; b2 = 1.348765 
  #random overlap
  # b1 = 0.00975   ; b2 = 1
    
  y = -100*(exp(b1*x^b2))^-1+100
  return(y)
}


cdfx <- cdfoverlap
#build models
# cdf2 <- subset(cdfx, crowns <= 700 & crowns >-10)
# 
# cdf2 <- subset(cdfx, crowns <= 295 & crowns >0)
# cv <- nls(tcov ~ -100*(exp(b1*crowns^b2))^-1+100, data = cdf2, 
#           start = list(b1 = 1, b2=1))
# cdf2 <- cdf2 |> mutate(tcov2 = predict(cv, crowns))
# ggplot(cdf2)+
#   geom_point(aes(x=crowns,y=tcov),color='red')+
#   geom_line(aes(x=crowns,y=tcov2))
# cv
# cvreverse <- nls(crowns ~  tcov*((1-b1*tcov/(tcov-100)))^b2, data = cdf2, 
#                  start = list(b1 = 1,b2 = 1))
# cdf2 <- cdf2 |> mutate(crowns2 = predict(cvreverse, tcov))
# ggplot(cdf2)+
#   geom_point(aes(x=crowns,y=tcov),color='red')+
#   geom_line(aes(x=crowns2,y=tcov))
# cvreverse

#test formulas
######
cdf2 <- subset(cdfavoid, crowns <= 295 & crowns >0)
cdf2 <- cdf2 |> mutate(k = tocov(crowns, avoid = T))
cv1 <- nls(crowns ~ k*((1-b1*k/(k-100)))^b2, 
           data = cdf2, start = list(b1 = 2.8950, b2 = 0.2834))#
cv2 <- nls(crowns ~ k*b3+b1*(log(100)-log(100-k))^b2, 
           data = cdf2, start = list(b1 = 1, b2 = 1, b3 = 1))#

cdf2 <- cdf2 |> mutate(kk1 = predict(cv1, k),
                       kk2 = predict(cv2, k))
ggplot(cdf2)+
  geom_point(aes(y=crowns,x=k),color='green')+
  geom_line(aes(y=crowns,x=crowns),color='yellow')+
  geom_line(aes(y=kk1,x=k),color='red')+
  geom_line(aes(y=kk2,x=k),color='blue')+
  scale_x_continuous(limits = c(0,100))

cdf2 <- cdf2 |> mutate(r1 = log(crowns) - log(kk1),
                       r2 = log(crowns) - log(kk2))
ggplot(cdf2)+
  geom_line(aes(y=r1,x=k),color='red')+
  geom_line(aes(y=r2,x=k),color='blue')


cdf3 <- subset(cdf2, k <= 90)
x0=cdf2$tcov
x1=cdf2$kk1
x2=cdf2$kk2

sum((x1-x0)^2)
sum((x2-x0)^2)

sum((log(x1)-log(x0))^2)
sum((log(x2)-log(x0))^2)




#get aggregate crown area
tocov <- function(kk, avoid=T){
  if(avoid){
    #avoid overlap
    b1=5.504e+01;b2=2.081e-02;b3=1.489e+02;b4=9.593e+03
    k = 2 * (100/(1 + exp(0 - kk/b1))^1 - 50) + b2 * kk*exp(-((kk - b3)^2/b4))
  }else{
    #random overlap
    b1=80.3178;b2=0.6273;b3=-124.8550;b4=23848.2169
    k = 2 * (100/(1 + exp(0 - kk/b1))^1 - 50) + b2 * kk*exp(-((kk - b3)^2/b4))
  }
  return(k)
}
#find crown area index
carea <- function(k, avoid=T){
  #k = canopy cover %
  if(avoid){
    #avoid overlap
    b1 = 23.2610; b2 = 1.1241; b3 = 0.9346
    kk = k * b3 + b1 * (log(100) - log(100 - k))^b2
  }else{
    #random overlap
    b1 = 88.3478; b2 = 1.0353; b3 = 0.2209
    kk = k * b3 + b1 * (log(100) - log(100 - k))^b2
  }
  #component area
  return(kk)
}

nstem <- function(k,cw, a=1, avoid=T){
  #k = canopy cover %
  #crown width m
  #a = area ha
  a0 = 10000*a
  #component area
  kk = carea(k=k, avoid=avoid)
  #relative crown area per unit area
  sa = (cw/2)^2*pi/a0
  #number of stems oer unit area
  st = round(kk/sa/100, 0)
  return(st)
}
#find crown width
findcw <- function(k, st, a=1, avoid=T){
  #k = canopy cover %
  #crown width m
  #a = area ha
  a0 = 10000*a
  #component area
  kk = carea(k=k, avoid=avoid)
  #relative crown area per unit area
  sa = a0*kk/100/st
  #number of stems oer unit area
  cw = round((sa/pi)^0.5*2,1)
  return(cw)
}
xxxxtocov <- function(kk, avoid=T){
  if(avoid){
    #avoid overlap
    b1 = 0.003639; b2 = 1.289286 
  }else{
    #random overlap
    b1 = 0.009502; b2 = 1.008376}
  k = -100*(exp(b1*kk^b2))^-1+100
  return(k)
}
#find crown area index
xxxxcarea <- function(k, avoid=T){
  #k = canopy cover %
  if(avoid){
    #avoid overlap
    b1 = 4.1225; b2 = 0.1507
  }else{
    #random overlap
    b1 = 2.8950; b2 = 0.2834}
  #component area
  kk = k*((1-b1*k/(k-100)))^b2
  return(kk)
}



nstem(k=50, cw=15, a=1)
findcw(k=50, st=32, a=1)

avoid=F
p=30
x <- carea(p,avoid = avoid)+carea(p,avoid = avoid)+carea(p,avoid = avoid)
x
y=tocov(x,avoid = avoid)
y
y=vegnasis::cover.agg(c(p,p,p))
y
cdf3 <- cdf

cp <- c('white',colorRampPalette(c('green','darkgreen'))(3))
stand$category <- cut(stand$shd,
                      breaks = c(-Inf,0.1, 0.2,Inf),
                      labels = c(0,0.1, 0.2),
                      
                      right = FALSE)
ggplot(stand, aes(x=xp,y=yp, color = category))+
  geom_point()+
  coord_fixed()+
  scale_color_manual(values =  cp)

ggplot(stand, aes(x=xp,y=yp, color = cd))+
  geom_point()+
  coord_fixed()+
  scale_color_gradient(low = 'green', high = 'darkgreen')



#totalcover

totcov <- subset(stand, xp >= 10 & xp <= 90 & yp >= 10 & yp <= 90) |> 
  mutate(shd =  ifelse(shd >0,1,0))
mean(totcov$shd)

# library(aqp)
# library(soilDB)
# colors <- c('10R 5/12', '10YR 9/12', '10Y 10/12', '10GY 8/12')
# weights <- c(60, 30, 10)
# 
# # combine into a data.frame and convert to sRGB + CIE LAB
# d <- cbind(
#   parseMunsell(colors, convertColors=FALSE),
#   parseMunsell(colors, return_triplets=TRUE, returnLAB=FALSE),
#   # pct=weights,
#   col=parseMunsell(colors, convertColors=TRUE)
# )
# 
# d <- d |> cbind(t(rgb2hsv(d$r*255, d$g*255, d$b*255)))
# d


# 
# col = data.frame(r=c(1,1,1,0,0,0,1), g=c(0,0.8,1,1,1,0,0), b=c(0,0,0,0,1,1,1))
# col = data.frame(h=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5)/6, s=c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5), v=c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5))
# col = c('brown', 'red', 'orange', 'yellow', 'green', 'blue','#4000FF',
#         '#8000FF','purple',
#         'violet','magenta','#800080')
# col2Munsell(col)
# col2rgb(col)|>t()
# rgb2hsv(col2rgb(col))|>t()*6
# col2Munsell(hsv(col$h, col$s, col$v))
