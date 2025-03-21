library(ggplot2)
library(vegnasis)

renumber <- function(x){
  n <- length(x)
  df <- data.frame(s=x, neworder = NA, ind = 1:n)
  dfsort <- df[order(df$s),]
  n2 <- length(x[!is.na(x)])
  if(n2>0){
  dfsort[!is.na(dfsort$s),]$neworder <- 1:n2
  s <- dfsort[order(dfsort$ind),]$neworder}else{s <- x}
  return(s)}



set.seed(2)
df <- data.frame(
  x=c(runif(50,0,50),0,50),
  y=c(rnorm(50,5,5),20,20),
  s = NA)
df <- data.frame(
  x=c(runif(50,0,50)),
  y=c(rnorm(50,5,15)),
  s = NA)



x=df$x
y=df$y
concavity = 1; curvy = TRUE; mag = 1

cavhull2 <- function(x,y, concavity = 0, curvy = FALSE, mag = 1, deep=FALSE){
  n=5 #number of segments to search between convex faces
  df <- data.frame(x=floor(x*1000)/1000,y=floor(y*1000)/1000) |> unique()
  
  #convex hull ----
  check <- TRUE #stopping rule
  #find starting point at bottom of plot
  miny <- min(df$y)
  minx <- min(df[df$y==miny,]$x)
  df <- mutate(df, s = ifelse(x==minx & y==miny, 1,NA))
  x1 <- df[df$s %in% 1,]$x
  y1 <- df[df$s %in% 1,]$y
  df <- df |> mutate(l1 = ((x-x1)^2+(y-y1)^2)^0.5,
                     a1=acos((x-x1)/l1)/2/pi*360,
                     a1=ifelse(y-y1 >=0,a1,-1*a1))
  amin = min(subset(df, !s %in% 1)$a1)
  df <- df |> mutate(s = ifelse(a1 == amin & is.na(s), 0, s))
  for(i in 1:nrow(df)){
    if(check){
      x0 <- df[df$s %in% (i-1),]$x
      y0 <- df[df$s %in% (i-1),]$y
      x1 <- df[df$s %in% i,]$x
      y1 <- df[df$s %in% i,]$y
      l0 = ((x1-x0)^2+(y1-y0)^2)^0.5
      a0 = acos((x1-x0)/l0)
      a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
      df <- df |> mutate(xr= x-x1,
                         yr= y-y1,
                         h=(xr^2+yr^2)^0.5,
                         a1=acos(yr/h),
                         a1=ifelse(xr >=0,a1,-1*a1),
                         a1= a1+a0,
                         xr = ifelse(h==0,0,h*sin(a1)),
                         yr = ifelse(h==0,0,h*cos(a1)),
                         xr = xr,
                         a1=asin(yr/h),
                         a1=ifelse(xr >=0,-a1,pi+a1))
      
      amin = min(subset(df, !s %in% c(i-1,i) )$a1)
      lmin = min(subset(df, !s %in% c(i-1,i) & a1 %in% amin)$l1)
      df <- df |> mutate(s = ifelse(s %in% 0,NA,s))
      check <- is.na(subset(df,a1 %in% amin & l1 %in% lmin)$s)
      df <- df |> mutate(s = ifelse(a1 == amin & is.na(s) & l1 %in% lmin, i+1, s))
    }
  }
  #concave hull first degree ----
  df <- df |> mutate(s1 = s, type = 'core')
  
  
  if(concavity > 0){
    for(k in 1:concavity){ #k=1
      smax <- max(df$s, na.rm = TRUE)
      #visit each convex hull boundary and rotate to a common reference
      for(i in 1:smax){#i=2
        i0 = ifelse(i == 1,smax,i-1)
        x0 <- df[df$s %in% (i0),]$x
        y0 <- df[df$s %in% (i0),]$y
        x1 <- df[df$s %in% i,]$x
        y1 <- df[df$s %in% i,]$y
        l0 <- ((x1-x0)^2+(y1-y0)^2)^0.5
        a0 = acos((x1-x0)/l0)
        a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
        dfr <- vegnasis::rotate(x=df$x, y=df$y, a=a0/2/pi*360, cx=x0, cy=y0)
        df <- df |> mutate(xr= dfr$x-x0,
                           yr= dfr$y-y0,
                           xs=NA,xa=NA,ys=NA,yl0=NA,ydiff=NA,microinc=NA)
        #use wave to select closest concave points
        en <- pmax(3,floor(pmin(n,l0/5)))*3
        #deep curve
        if(deep & k == concavity){
          wave0 <- data.frame(x=(0:(en+1))/(en+1))
          wave0 <- wave0 |> mutate(a=x*2*pi,y=(cos(a)^1-1)/2*mag)
          wave <- data.frame(x=NA, y=NA, s=NA,l1=NA,a1=NA,xr=wave0$x*l0,
                             yr=wave0$y*l0,h=NA,s1=NA, type='cave',
                             xs=NA,xa=NA,ys=NA,yl0=NA,ydiff=NA,microinc=NA)
          wavr <- vegnasis::rotate(x=wave$xr, y=wave$yr, a=-a0/2/pi*360, cx=0,cy=0)
          wave <- wave |> mutate(x=wavr$x+x0,y=wavr$y+y0) |> subset(!yr >=0)
          df <- df |> rbind(rbind(wave))
        }
        df <- df |> mutate(xs = xr/l0, xa = xs*2*pi, ys = (cos(xa)^1-1)/2*mag,
                           yl0 = (yr/l0), ydiff = yl0-ys)
        curmax <- max(subset(df, xs > 0 & xs < 1)$ydiff)
        curcur <- subset(df, xs >= 0 & xs <=1  & ydiff == curmax)$ys
        currat <- ifelse(curmax == 0, 1, 1-curmax/abs(curcur))
        currat <- ifelse(currat > 1,1,currat)
        df <- df |> mutate(ys = ys*currat, ydiff = yl0-ys)
        curmax <- max(subset(df, xs > 0 & xs < 1)$ydiff)
        curcur <- subset(df, xs >  0 & xs < 1  & ydiff == curmax)$ys
        currat2 <- ifelse(curmax == 0, 1, 1-curmax/abs(curcur))
        currat2 <- ifelse(currat2 > 1,1,currat2)
        df <- df |> mutate(ys = ys*currat2, ydiff = yl0-ys)
        df$microinc <- NA
        #introduce wave
        if(curvy & k == concavity){
        wave0 <- data.frame(x=(0:(en+1))/(en+1))
        wave0 <- wave0 |> mutate(a=x*2*pi,y=(cos(a)^1-1)/2*mag)
        wave <- data.frame(x=NA, y=NA, s=NA,l1=NA,a1=NA,xr=wave0$x*l0,
        yr=wave0$y*l0*currat*currat2,h=NA,s1=NA, type='wave',
        xs=NA,xa=NA,ys=NA,yl0=NA,ydiff=NA,microinc=NA)
        wavr <- vegnasis::rotate(x=wave$xr, y=wave$yr, a=-a0/2/pi*360, cx=0,cy=0)
        wave <- wave |> mutate(x=wavr$x+x0,y=wavr$y+y0) |> subset(!yr >=0)
        df <- df |> rbind(rbind(wave))
        }
        pickthispoint <- min(abs(subset(df, xs >= 0 & xs <=1)$ydiff))
        df <- df |> mutate(microinc = ifelse((xs >  0 & xs < 1  & round(abs(ydiff),10) %in% round(pickthispoint,10) | type %in% 'wave') & is.na(s1),xr,NA))
       
        df$microinc <- renumber(df$microinc)
        
        df <- df |> mutate(s1 = ifelse(!is.na(microinc) & is.na(s), i0+microinc/1000,s1))
       df <- df |> subset(type %in% 'core' | !is.na(s1)) 
      }
            df <- df |> mutate(s = renumber(s1), s1 = s)
    }
  }
  df <- subset(df, !is.na(s), select=c(x,y,s)) |> arrange(s)
  return(df)}

d1 <- cavhull2(df$x,df$y, concavity = 0)
d2 <- cavhull2(df$x,df$y, concavity = 1, mag = 1)
d3 <- cavhull2(df$x,df$y, concavity = 2, curvy = F, mag = 1, deep=F)

ggplot()+
  geom_polygon(data=d1,aes(x=x,y=y), color='blue', fill='blue')+
  geom_polygon(data=d2,aes(x=x,y=y), color='red', fill='red')+
  geom_polygon(data=d3,aes(x=x,y=y), color='green', fill='green')+
  geom_point(data=df,aes(x=x,y=y))+
  coord_fixed()




##########################3 circles
#l0 <- 39.65157
pmc <- .3

h1 <- pmc*l0/2#center circle
h2 <- ((l0/2)^2-h1^2)/(2*h1)
cir <- data.frame(a=(1:360)/360*2*pi)
cir <- cir |> mutate(x=cos(a),y=sin(a))
a2 <- asin(h2/(h1+h2))
a2/2/pi*360
180-a2/2/pi*360-90
ccir <- cir |> subset(y < 0 & 2*pi-a >= a2 & a-pi >= a2) |> mutate(az=a-pi,x=x*h1,y=y*h1)
lcir <- cir |> subset(y > 0 & x > 0 & a >= a2) |> mutate(x=x*h2-l0/2,y=y*h2-h2)
rcir <- cir |> subset(y > 0 & x < 0 & pi-a >= a2) |> mutate(x=x*h2+l0/2,y=y*h2-h2)

ggplot()+
  geom_path(data=ccir,aes(x=x,y=y))+
  geom_path(data=lcir,aes(x=x,y=y))+
  geom_path(data=rcir,aes(x=x,y=y))+
  coord_fixed()

########################## sine wave
#l0 <- 39.65157
pmc <- .3

npts <- 8
wave <- data.frame(x=(0:npts)/npts)

wave <- wave |> mutate(a=x*2*pi,y=cos(a))
wave0 <- wave |> mutate(y=(y+1)/2)
wave1 <- wave0
wave2 <- wave |> mutate(y=-1*(y-1)/2)
wave1 <- wave1 |> mutate(y=y^2)
wave2 <- wave2 |> mutate(y=-1*y^2+1)



ggplot()+
  geom_path(data=wave0,aes(x=x,y=y))+
  geom_path(data=wave1,aes(x=x,y=y), color='red')+
  geom_path(data=wave2,aes(x=x,y=y), color='blue')+
  coord_fixed()



ggplot()+
  geom_point(data=subset(df2),aes(x=xs,y=yl0))+
  geom_path(data=arrange(subset(df, xs >= 0 & xs <=1), xs),aes(x=xs,y=ys), color='blue')
coord_fixed()