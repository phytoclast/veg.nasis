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
  

