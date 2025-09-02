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
cw = 7
scf = cw/9 
stand <- make_hex_stand(plength/100/scf,1*scf)
# nudgx <- runif(nrow(stand), min=-0.25, max=0.25)
# nudgy <- runif(nrow(stand), min=-0.25, max=0.25)
# stand <- stand |> mutate(xp=xp+nudgx, yp=yp+nudgy)
stand <- stand |> mutate(dbh = 0, cwd = 0, shd=0, d=NA, cd = 100)
for(i in 1:1000){#i=1
empty <- stand[stand$dbh == 0,]
empty <- empty[sample(1:nrow(empty),size=1, replace = FALSE, prob = empty$wt),]$stumpid
stand <- stand |> mutate(dbh = ifelse(stumpid %in% empty, 30,dbh),
                           cwd = ifelse(stumpid %in% empty, 9,cwd))
newtree <- subset(stand, stumpid %in% empty)
newtree <- newtree |> mutate(cwd = pmax(cwd/15,pmin(cwd, cd*15)))
stand <- stand |> mutate(d = ((newtree$yp-yp)^2+(newtree$xp-xp)^2)^0.5,
                           shd = ifelse(d < newtree$cwd/2, shd+0.1,shd),
                           cd = pmin(d-newtree$cwd/2, cd, na.rm = TRUE),
                           wt = ifelse(shd <= 0, 1,ifelse(shd <= 0.1, 0.5, 0.25)))
}


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
