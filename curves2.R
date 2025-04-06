library(ggplot2)
library(vegnasis)
df=readRDS('C:/scripts/veg.nasis/df.rds')
x=c(1,2,3,4,5, 4,3,2,1, 1.5,2.5,3.5,4.5, 1.6);
y=c(1,1,1,1,1, 2,3,2,1, 1.29,2.75,3.2,1.1, 1.27)
concavity = 1; curvy = T; mag = 1; maxdepth=0.7; minspan=0
#new inner border maker
df=data.frame(x=x,y=y)


x=df$x;y=df$y
# cavhull2 <- function(x,y, concavity = 0, curvy = FALSE, maxdepth=NA, minspan=0, mag = 1){
n=5 #number of segments to search between convex faces
concavity = 1; curvy = T; mag = 1; maxdepth=0.7; minspan=0
df <- data.frame(x=floor(x*100)/100,y=floor(y*100)/100) |> unique()

#convex hull ----
check <- TRUE #stopping rule
#find starting point at bottom of plot
miny <- min(df$y)
minx <- min(df[df$y==miny,]$x)
df <- mutate(df, s = ifelse(x==minx & y==miny, 1,NA))
x1 <- df[df$s %in% 1,]$x
y1 <- df[df$s %in% 1,]$y
df <- df |> mutate(l1 = ((x-x1)^2+(y-y1)^2)^0.5,
                   a1=acos((x-x1)/l1),
                   a1=ifelse(y-y1 >=0,a1,-1*a1))
amin = min(subset(df, !s %in% 1)$a1)
lmin = min(subset(df, !s %in% 1 & a1 %in% amin)$l1)
df <- df |> mutate(s = ifelse(a1 == amin & l1 == lmin & is.na(s), 0, s))
for(i in 1:nrow(df)){#nrow(df)
  if(check){#i=2
    x0 <- df[df$s %in% (i-1),]$x
    y0 <- df[df$s %in% (i-1),]$y
    x1 <- df[df$s %in% i,]$x
    y1 <- df[df$s %in% i,]$y
    l0 = ((x1-x0)^2+(y1-y0)^2)^0.5
    a0 = round(acos((x1-x0)/l0),9)
    a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
    df <- df |> mutate(xr= x-x1,
                       yr= y-y1,
                       l1=round((xr^2+yr^2)^0.5,9),
                       a1=round(acos(yr/l1),9),
                       a1=ifelse(xr >=0,a1,-1*a1),
                       a1= a1+a0,
                       xr = round(ifelse(l1==0,0,l1*sin(a1)),9),
                       yr = round(ifelse(l1==0,0,l1*cos(a1)),9),
                       a1=asin(round(yr/l1,9)),
                       a1=ifelse(xr >=0,-a1,pi+a1))
    
    amin = min(subset(df, !s %in% c(i-1,i) )$a1)
    lmin = min(subset(df, !s %in% c(i-1,i) & a1 %in% amin)$l1)
    df <- df |> mutate(s = ifelse(s %in% 0,NA,s))
    check <- is.na(subset(df,a1 %in% amin & l1 %in% lmin)$s)
    df <- df |> mutate(s = ifelse(a1 == amin & is.na(s) & l1 %in% lmin, i+1, s))
  }
}

df <- df |> mutate(s1 = s, type = 'core')
smax <- max(df$s, na.rm = TRUE)
#visit each convex hull boundary and rotate to a common reference
# for(i in 1:smax){#
  i=6
  i0 = ifelse(i == 1,smax,i-1)
  x0 <- df[df$s %in% (i0),]$x
  y0 <- df[df$s %in% (i0),]$y
  x1 <- df[df$s %in% i,]$x
  y1 <- df[df$s %in% i,]$y
  l0 <- ((x1-x0)^2+(y1-y0)^2)^0.5
  # if(l0 > minspan){
    a0 = acos((x1-x0)/l0)
    a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
    dfr <- vegnasis::rotate(x=df$x, y=df$y, a=a0/2/pi*360, cx=x0, cy=y0)
    df <- df |> mutate(xr= dfr$x-x0,
                       yr= dfr$y-y0)
    #-----------
    
    df <- df |> mutate(xs = round(xr/l0,6))
    n = 20
    
    qymax <- matrix(NA, nrow = n, ncol = 5)
    
    for(q in 1:n){#q=2
      px0 = (q-1)/n*l0
      px1 = (q)/n*l0
      qymax[q,1] <- max(c(-1,subset(df, xr >0 & xr >= px0 & xr < px1)$yr))
      qx <- subset(df, xr >0 & xr >= px0 & xr < px1 & yr %in% qymax[q,1])$xr[1]
      qymax[q,2] <- ifelse(is.na(qx), (px0+px1)/2, qx)
     }
    for(q in 1:n){#q=1
      if(q==1){
        qymax[q,3] <- 0
      }else{
        qymax[q,3] <- pmax(qymax[q,1],qymax[q-1,3]-0.02*l0)
      }
    }
    for(q in n:1){#q=20
      if(q==n){
        qymax[q,4] <- 0
      }else{
        qymax[q,4] <- pmax(qymax[q,1],qymax[q+1,4]-0.02*l0)
      }
    }
    qymax[,5] <- pmax(qymax[,3],qymax[,4])
    
    qymax <-as.data.frame(qymax)
    
    
    
    
    
    
    
  # }}
    # df2 <- subset(df, !is.na(s), select=c(x,y,s)) |> arrange(s)
    df2 <- subset(df, !is.na(s)) |> arrange(s)
  
  ggplot()+
    geom_polygon(data=df2,aes(x=xr,y=yr), color='lightblue', fill='lightblue')+
    geom_point(data=qymax,aes(x=V2,y=V5), color='green', fill='green')+
    geom_point(data=df,aes(x=xr,y=yr), color='red', fill='red')+
    geom_point(data=df2,aes(x=xr,y=yr), color='black', fill='black')+
    coord_fixed()
  
  ggplot()+
    geom_polygon(data=df2,aes(x=x,y=y), color='lightblue', fill='lightblue')+
    geom_point(data=df,aes(x=x,y=y), color='red', fill='red')+
    geom_point(data=df2,aes(x=x,y=y), color='black', fill='black')+
    coord_fixed()
  
