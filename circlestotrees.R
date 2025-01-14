library(ggplot2)
library(dplyr)
angles <- c(0:599)/600*2*3.141592
#coniferoverhead

x = cos(angles)/2; y=sin(angles)/2
ht=3
nbr = 8*6/4/ht^0.75
tree <- data.frame(shape='shape', x=round(x,5),y=round(y,5))
tree <- tree |> mutate(y = (y+0.2)/1.2*ht)
tree2 <- tree |> mutate(shape='shape2')
tree2 <- tree2 |> mutate(angles=angles/2/3.141592*360, 
                         top = max(y), d = round(((((y-top)^2+x^2)^0.5)*nbr),3), 
                         d1 = floor(d), yy = (y>=0)) |> group_by(d1, yy) |> 
  mutate(dmin = min(d)) |> ungroup() |> mutate(p = (dmin==d))  |> subset(y >= 0)
top <- tree2$top[1]
for(i in 1:15){
  treex0 <- tree |> mutate(circle=i, x=cos(angles)/nbr*i, y=sin(angles)/nbr*i+top)
  if(i==1){treex=treex0}else{treex<-rbind(treex,treex0)}}
tree2 <- subset(tree2, p)

ggplot()+
  geom_polygon(data=tree, aes(x=x, y=y), fill='#FFFFFF00')+
  geom_point(data=tree, aes(x=x, y=y), color='gray')+
  geom_point(data=tree2, aes(x=x, y=y), color='red')+
  # geom_polygon(data=treex, aes(x=x, y=y, group=circle), color='yellow', fill='#FFFFFF00')+
  coord_fixed(ratio = 1)

######
library(ggplot2)
library(dplyr)



df <- data.frame(i = 0:10)
df <- df |> mutate(s = 1-i/10,  x = s^0.7, y = 1-(s^3))
df <- df |> mutate(i0 = i+0, i1 = i+0.3, i2 = i+0.7, 
                   x1 = (s+1)*0.02, y1 = (1-s)*0.50)


df.tips <- data.frame(i = df$i1, x = df$x, y = df$y)
df.base1 <- data.frame(i = df$i0, x = df$x1, y = df$y1-0.01)
df.base2 <- data.frame(i = df$i2, x = df$x1, y = df$y1+0.01)

tree <- rbind(df.base1,df.tips,df.base2) |> arrange(i)
tree <- tree[-nrow(tree),]
tree2 <- tree |> mutate(i= -i+21, x=-x)
tree2 <- rbind(tree, tree2) |> arrange(i)

ggplot()+
  geom_point(data=df, aes(x=x, y=y), color='red')+
  geom_point(data=df, aes(x=x1, y=y1), color='red')+
  geom_polygon(data=tree2, aes(x=x, y=y), color='red')+
  coord_fixed()

######
library(ggplot2)
library(dplyr)


makeStem <- function(lth, wth, tip=0.01, inc=10){
  
  
  s <- (0:inc)/inc
  i <- c(1:(inc+1),(1:(inc+1))+inc+1)
  y <- s*lth; y <- c(y,y[(inc+1):1])
  x <- tip/2*s+wth/2*(s*-1+1); x <- c(-1*x,x[(inc+1):1])
  side <- ifelse(x >= 0, 'R','L')
  df <- data.frame(x = x,
                   y = y,
                   i = i,
                   type = 'mid',
                   side=side,
                   center=0,
                   width=abs(x*2)
  )
  df[c(1,nrow(df)),]$type <- 'base'
  df[c(nrow(df)/2,nrow(df)/2+1),]$type <- 'tip'
  return(df)
}


angle=30
bht=2
lth=4
wth = 1
sc = 0.6
stem <-  makeStem(lth,wth,0.1*wth,20)
branch <-  makeStem(0.2,0.1,0.01,10)
stem <-  stem |> mutate(x = x + -0.2*cos(y/lth*2*pi), center = center + -0.2*cos(y/lth*2*pi))
branch <-  stem |> mutate(x=x*sc,y=y*sc,center=center*sc, width = width*sc)
branch2 <-  stem |> mutate(x=x*sc*0.7,y=y*sc*0.7,center=center*sc*0.7)
stem <- attachBranch(stem, branch, -50, 2)
branch <- branch2


attachBranch <- function(stem, branch, angle, bht){
  #establish permanent columns to end up with
  original <- colnames(stem)
  
  #get information about stem length
  stemax <- max(subset(stem, type %in% c('base','mid','tip'))$y)
  stemin <- min(subset(stem, type %in% c('base','mid','tip'))$y)
  
  #ensure that branch never rises higher than stem
  bht <- ifelse(bht > stemax,stemax,bht)
  
  #ensure that branch is no thicker than stem
  bhtlower <- max(stem[bht >= stem$y,]$y) 
  bhtupper <- min(stem[bht <= stem$y,]$y)
  swd <- (mean(subset(stem, y == bhtlower)$width)*1/(abs(bht-bhtlower)+0.01)+
            mean(subset(stem, y == bhtupper)$width)*1/(abs(bht-bhtupper)+0.01))/(1/(abs(bht-bhtlower)+0.01)+1/(abs(bht-bhtupper)+0.01))
  bbase <- subset(branch, type %in% 'base')
  bwd <- max(bbase$x)-min(bbase$x)
  if(bwd > swd) {branch <- branch |> mutate(x = x*swd/bwd)}
  
  #convert angle to radians
  angle = angle/360*2*pi
  
  #set which side of stem branch will be fitted
  xside <- ifelse(angle>0,'R','L')
  xsine <- ifelse(angle>0,1,-1)
  
  #rotate branch to correct angle
  branch <- branch |> mutate(h=(x^2+y^2)^0.5,a = acos(y/h),a=ifelse(x < 0,-1*a,a),
                             x = h*sin(a+angle), y = h*cos(a+angle))
  
  #lift branch to correct height
  branch <- branch |> mutate(y = y+bht)
  
  #recheck to see if branch thickness pushes branch too high
  bdif <- stemax - max(branch[branch$type %in% 'base',]$y)
  if(bdif < 0){branch <- branch |> mutate(y = y+bdif)}
  
  #determine which stem vertices straddle the branch
  xd <- subset(stem, side %in% xside & type %in% c('tip','mid','base','bbase'))
  ylower <- max(subset(xd, y < min(branch[branch$type %in% 'base',]$y))$y)
  yupper <- min(subset(xd, y >= max(branch[branch$type %in% 'base',]$y))$y)
  xupper <- subset(xd, y == yupper)$x
  xlower <- subset(xd, y == ylower)$x
  
  #shift branch to conform with twisted stem center position
  cshift <- mean(subset(xd, y == yupper | y == ylower)$center)
  branch <- branch |> mutate(x = x+cshift)
  
  #identify which branch vertices are inside stem
  branch <- branch |> mutate(inside = (x-((y-ylower)/(yupper-ylower)*(xupper-xlower)+xlower))*xsine)
  
  #identify vertices which straddle inside and outside of stem to find points of intersection
  internal <- branch |> mutate(near = inside^2, isinside = ifelse(inside < 0 | type %in% 'base', 'no','yes')) |> group_by(isinside, side) |> mutate(minnear = min(near)) |> ungroup() |> subset(near == minnear)
  
  #approximate location of branch stem intersection to establish new branch base
  newbase <- internal |> group_by(side) |> mutate(amt = 1/((inside - 0)/(max(inside)-min(inside)+0.0001))^2) |> 
    summarise(x= sum(amt*x)/sum(amt), y= sum(amt*y)/sum(amt), i= mean(i), type='base', center=0, side=xside, width=bwd) 
  
  
  #determine if stem if branch too close to top
  if(stemax - bht < 0.05*(stemax-stemin)){
    #identify where to insert new numbering sequence to maintain correct vertex order
    xdi <- mean(subset(stem, type %in% 'tip')$i)
    #remove tip of stem
    steminternal <- stem |> subset(!type %in% 'tip' & !(y > (stemax - 0.05*(stemax-stemin)) & type %in% 'mid'))
    #assemble branch with new base, omitting internal vertices
    branchinternal <- branch |> subset(select=original) |> rbind(newbase[,original]) |> mutate(i = xdi + i/10000, type=paste0('b',type), inside = NULL, h=NULL,a=NULL)  |> arrange(i) 
  }else{
    #remove stem vertices that may be covered by new branch
    steminternal <- stem |> subset(!(x >= min(internal$x) & x <= max(internal$x) & y >= min(internal$y) & y <= max(internal$y))) |> subset(select=original)
    #identify stem vertices near branch base
    ylower2 <- max(subset(xd, y < min(newbase$y))$y)
    yupper2 <- min(subset(xd, y > max(newbase$y))$y)
    #identify where to insert new numbering sequence to maintain correct vertex order
    xdi <- mean(subset(xd, y %in% c(yupper2, ylower2))$i)
    #assemble branch with new base, omitting internal vertices
    branchinternal <- branch |> subset(inside >= 0) |> subset(select=original) |> rbind(newbase[,original]) |> mutate(i = xdi + i/10000, type=paste0('b',type), inside = NULL, h=NULL,a=NULL)  |> arrange(i) 
    
  }
  
  
  #append branch to stem with correct vertex order
  stemnew <- rbind(branchinternal,steminternal) |> arrange(i) 
  
  #renumber vertices
  stemnew <- mutate(stemnew, i=(1:nrow(stemnew)))
  
  return(stemnew)}

ggplot()+
  geom_polygon(data=stemnew, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=stemnew, aes(x=x, y=y), color='red')+
  geom_polygon(data=branch, aes(x=x, y=y), color='blue',fill='#00009950')+
  geom_point(data=branch, aes(x=x, y=y), color='blue')+
  # geom_polygon(data=branchinternal, aes(x=x, y=y), color='green',fill='#00990050')+
  # geom_point(data=branchinternal, aes(x=x, y=y), color='green')+
  coord_fixed()

skewStem <- function(stem, amp=0.2, phase=0, waves=1){
  maxstem <- max(stem$y)
  minstem <- min(stem$y)
  lth <- maxstem - minstem
  stem <-  stem |> mutate(x = x + amp*cos((y/lth+phase)*2*pi*waves), 
                          center = center + amp*cos((y/lth+phase)*2*pi*waves))
  return(stem)
}

stem <-  makeStem(lth,wth,0.2*wth,20)
stem <-  skewStem(stem,-0.1,0.25,3)


ggplot()+
  geom_polygon(data=stem, aes(x=x, y=y), color='brown',fill='#99500050')+
  coord_fixed()
#------
#fractal tree with bend --- 
n=3
lth=4
wth = 0.5
sc = 0.6
stem <-  makeStem(lth,wth,0.1*wth,20)
# stem <-  stem |> mutate(x = x + -0.2*cos(y/lth*2*pi), center = center + -0.2*cos(y/lth*2*pi))
stem <-  skewStem(stem,-0.2,.5,3)
branch <-  stem |> mutate(x=x*sc,y=y*sc,center=center*sc)
branch2 <-  stem |> mutate(x=x*sc*0.7,y=y*sc*0.7,center=center*sc*0.7)
branch3 <-  stem |> mutate(x=x*sc*0.3,y=y*sc*0.3,center=center*sc*0.3)
for(i in 1:n){
  tree <- attachBranch(stem, branch, -50, 2)
  branch <- branch2
  tree <- attachBranch(tree, branch, 50, 3)
  branch <- branch3
  #stem=tree
  tree <- attachBranch(tree, branch, -5, 4)
  branch <- tree |> mutate(x=x*sc,y=y*sc,center=center*sc)
  branch2 <- tree |> mutate(x=x*sc*.7,y=y*sc*.7,center=center*sc*.7)
  branch3 <- tree |> mutate(x=x*sc*.3,y=y*sc*.3,center=center*sc*.3)
}

crown <- tree |> subset(grepl('tip',type) & !type %in% 'tip')
ggplot()+
  geom_polygon(data=tree, aes(x=x, y=y), color='brown',fill='#99500050')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=crown, aes(x=x, y=y), color='green',fill='#00990050')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()
#------


makeCrowShape <- function(n=9, h, w, shape=c('pyramid','dome','round','column'))
  shape=c('pyramid','dome','round','column')
h=2
w=5
wd=w/2
n=6

r <- c('pyramid','dome','round','column') %in% shape |> as.numeric(r)

s = (0:n)/n
i = 1:(n+1)
d = 1-wd/h
sc <- (s-d)*1/(1-d)
p=0
opposite = TRUE
shapes <- data.frame(i=i,s=s)
shapes <- shapes |> mutate(sc = (s-d)*1/(1-d),
                           a1=s*2*pi/4,
                           a2=(s*2-1)*2*pi/4,
                           a3=sc*2*pi/4)

shapes <- shapes |> mutate(x0=-1*s+1,y0=s,#pyramid
                           x1=cos(a1),y1=sin(a1),#dome
                           x2=cos(a2),y2=sin(a2),#round 
                           x3=cos(a3),y3=sin(a3),#column
                           x3=ifelse(s > d,x3,1),y3=ifelse(s > d,(y3*(1-d)+d),y0))
shapes <- shapes |> mutate(tx = (r[1]*x0+r[2]*x1+r[3]*x2+r[4]*x3)/sum(r),
                           ty = (r[1]*y0+r[2]*y1+r[3]*y2+r[4]*y3)/sum(r))
#normalize shape
shapes <- shapes |> mutate(ty=h*(ty-min(ty))/(max(ty)-min(ty)), tx=wd*(tx-min(tx))/(max(tx)-min(tx)))
#branch base
shapes <- shapes |> mutate(bx = 0, by=(h-w)*s) |> subset(select = c(i,s,tx,ty,bx,by))
shapes <- shapes |> mutate(a = 360/(2*pi)*acos((ty-by)/((ty-by)^2+(tx-bx)^2)^0.5))
if(opposite){
  shapes2 <- shapes |> mutate(a = a*-1, i=i+0.5) |> subset(!a %in% 0)
  shapes <- shapes |> rbind(shapes2) |> arrange(i) 
  shapes <- shapes |> mutate(i = 1:nrow(shapes))
}else{
  shapes <- shapes |> mutate(a = ifelse(i/2 == floor(i/2), a*-1,a))
}

ggplot()+
  # geom_polygon(data=circle, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=shapes, aes(x=bx, y=by), color='red')+
  # geom_point(data=shapes, aes(x=x0, y=y0), color='blue')+
  geom_point(data=shapes, aes(x=tx, y=ty), color='green')+
  geom_path(data=shapes, aes(x=bx, y=by), color='red')+
  # geom_path(data=shapes, aes(x=x0, y=y0), color='blue')+
  geom_path(data=shapes, aes(x=tx, y=ty), color='green')+
  # geom_point(data=shapes2, aes(x=x, y=y), color='purple')+
  # geom_point(data=circle2, aes(x=x, y=y), color='gold')+
  coord_fixed()


makeCrowShape <- function(ht.max=5, ht.min=1, crwd=2, dbh, crshape=c('pyramid','dome','round','column'), n=5, bu=0.8, bl=0, opposite = FALSE){
  
  h <- ht.max - ht.min
  wd <- crwd/2
  
  s <- (0:n)/n
  i <- 1:(n+1)
  d <- 1-wd/h
  sc <- (s-d)*1/(1-d)
  #create index sequence and proportion
  shapes <- data.frame(i=i,s=s)
  #create angles needed for round crowns
  shapes <- shapes |> mutate(sc = (s-d)*1/(1-d),
                             a1=s*2*pi/4,
                             a2=(s*2-1)*2*pi/4,
                             a3=sc*2*pi/4)
  #create choices of xy coordinates
  shapes <- shapes |> mutate(x0=-1*s+1,y0=s,#pyramid
                             x1=cos(a1),y1=sin(a1),#dome
                             x2=cos(a2),y2=sin(a2),#round 
                             x3=cos(a3),y3=sin(a3),#column
                             x3=ifelse(s > d,x3,1),y3=ifelse(s > d,(y3*(1-d)+d),y0))
  #weights to determine shape output as branch tip coordinates
  r <- c('pyramid','dome','round','column') %in% crshape |> as.numeric()
  shapes <- shapes |> mutate(tx = (r[1]*x0+r[2]*x1+r[3]*x2+r[4]*x3)/sum(r),
                             ty = (r[1]*y0+r[2]*y1+r[3]*y2+r[4]*y3)/sum(r))
  #normalize shape
  shapes <- shapes |> mutate(ty=h*(ty-min(ty))/(max(ty)-min(ty)), tx=wd*(tx-min(tx))/(max(tx)-min(tx)))
  #branch base
  shapes <- shapes |> mutate(bx = 0, by=(bu*h-bl*h)*s+bl*h) |> subset(select = c(i,s,tx,ty,bx,by))
  #lift branches to crown base
  shapes <- shapes |> mutate(by=by+ht.min, ty=ty+ht.min)
  #angle of branch
  shapes <- shapes |> mutate(l = pmax(dbh/2,((ty-by)^2+(tx-bx)^2)^0.5), a = 360/(2*pi)*acos((ty-by)/l))
  #branch diameter
  shapes <- shapes |> mutate(d = dbh*pmin(l*2,(ht.max-by))/ht.max+0.01)
  if(opposite){
    shapes2 <- shapes |> mutate(a = a*-1, tx = tx*-1, i=i+0.5) |> subset(!a %in% 0)
    shapes <- shapes |> rbind(shapes2) |> arrange(i) 
    shapes <- shapes |> mutate(i = 1:nrow(shapes))
  }else{
    shapes <- shapes |> mutate(o = ifelse(i/2 == floor(i/2), -1,1), a = a*o, tx = tx*o, o = NULL)
  }
  
  return(shapes)
}

shapes <- makeCrowShape(ht.max=10, ht.min=4, crwd=6, dbh=0.5, crshape=c('pyramid'), n=10, bu=0.8, bl=-0.2, opposite = F)

stem <-  makeStem(10,1,0.1,10)
for(i in 1:nrow(shapes)){#i=1
  branch <- makeStem(shapes$l[i], shapes$d[i],0.01,10)
  branch <- skewStem(branch, amp=ifelse(shapes$a[i] >= 0,-0.2,0.2), phase=0, waves=1.2)
  stem <- attachBranch(stem, branch, shapes$a[i], shapes$by[i])
}
crown <- stem |> subset(grepl('tip',type))


ggplot()+
  geom_polygon(data=stem, aes(x=x, y=y), color='brown',fill='#99500050')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=crown, aes(x=x, y=y), color='green',fill='#00990050')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()



ggplot()+
  # geom_polygon(data=circle, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=shapes, aes(x=bx, y=by), color='red')+
  # geom_point(data=shapes, aes(x=x0, y=y0), color='blue')+
  geom_point(data=shapes, aes(x=tx, y=ty), color='green')+
  geom_path(data=crown, aes(x=bx, y=by), color='red')+
  # geom_path(data=shapes, aes(x=x0, y=y0), color='blue')+
  geom_path(data=crown, aes(x=tx, y=ty), color='green')+
  # geom_point(data=shapes2, aes(x=x, y=y), color='purple')+
  # geom_point(data=circle2, aes(x=x, y=y), color='gold')+
  coord_fixed()








stem <-  makeStem(20,1,0.3*wth,20)
for(i in 1:6){#i=2
  lth = dfbranches$crownht[i]
  bht <- dfbranches$baseht[i]
  angle = dfbranches$angle[i] * ifelse(i/2==floor(i/2),-1,1)
  branch <- makeStem(lth,lth/8,lth/13)
  stem <- attachBranch(stem,branch,angle,bht)
  
}


x <- stem %>%
  st_as_sf(coords = c("x", "y")) |>
  summarise(i = list(i), type = list(type),
            side = list(side), center = list(center), 
            geometry = st_combine(geometry)) |> st_cast("POLYGON")


x <- data.frame(x=st_coordinates(x)[,1],
                y=st_coordinates(x)[,2],
                i=unlist(x$i),
                type=unlist(x$type),
                side=unlist(x$side),
                center=unlist(x$center))

poly_sf <- function(x){
  require(sf)
  x <- x |> st_as_sf(coords = c("x", "y")) |>
    summarise(i = list(i), type = list(type),
              side = list(side), center = list(center), 
              geometry = st_combine(geometry)) |> st_cast("POLYGON")
  return(x)
}
sf_poly <- function(x){
  require(sf)
  x <- data.frame(x=st_coordinates(x)[1:(nrow(x)-1),1],
                  y=st_coordinates(x)[1:(nrow(x)-1),2],
                  i=unlist(x$i),
                  type=unlist(x$type),
                  side=unlist(x$side),
                  center=unlist(x$center))
  return(x)
}

x <- poly_sf(stem)

stem2 <- sf_poly(x)
angle=50
bht=3
lth=4
wth = 0.5
sc = 0.6

stem <-  makeStem(lth,wth,0.1*wth,20)
stem <-  stem |> mutate(x = x + -0.2*cos(y/lth*2*pi), center = center + -0.2*cos(y/lth*2*pi))
branch <-  stem |> mutate(x=x*sc,y=y*sc,center=center*sc)













tree2 <- data.frame(x=st_coordinates(polygon)[,1],y=st_coordinates(polygon)[,2])
ggplot()+
  geom_polygon(data=tree2, aes(x=x, y=y), color='green',fill='#00990050')+
  coord_fixed()