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


makeStem <- function(lth, wth, tip, inc=10){
  
  
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


angle=-30
bht=4
lth=4
wth = 0.5
sc = 0.6
stem <-  makeStem(lth,wth,0.1*wth,20)
stem <-  stem |> mutate(x = x + -0.2*cos(y/lth*2*pi), center = center + -0.2*cos(y/lth*2*pi))
branch <-  stem |> mutate(x=x*sc,y=y*sc,center=center*sc, width = width*sc)
branch2 <-  stem |> mutate(x=x*sc*0.7,y=y*sc*0.7,center=center*sc*0.7)
stem <- attachBranch(stem, branch, -50, 2)
branch <- branch2
  

attachBranch <- function(stem, branch, angle, bht){
  #establish permanent columns to end up with
  original <- colnames(stem)
  #ensure that branch never rises higher than stem
  bht <- ifelse(bht > max(stem$y),max(stem$y),bht)
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
  bdif <- max(subset(stem, type %in% c('base','mid','tip'))$y) - max(branch[branch$type %in% 'base',]$y)
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
  
  #? bx <- (bht-ylower)/(yupper-ylower)*(xupper-xlower)+xlower
  #identify which branch vertices are inside stem
  branch <- branch |> mutate(inside = (x-((y-ylower)/(yupper-ylower)*(xupper-xlower)+xlower))*xsine)
  #identify vertices which straddle inside and outside of stem to find points of intersection
  internal <- branch |> mutate(near = inside^2, isinside = ifelse(inside < 0 | type %in% 'base', 'no','yes')) |> group_by(isinside, side) |> mutate(minnear = min(near)) |> ungroup() |> subset(near == minnear)
  #approximate location of branch stem intersection to establish new branch base
  newbase <- internal |> group_by(side) |> mutate(amt = 1/((inside - 0)/(max(inside)-min(inside)))^2) |> 
    summarise(x= sum(amt*x)/sum(amt), y= sum(amt*y)/sum(amt), i= mean(i), type='base', center=0, side=xside, width=bwd) 
  #remove stem vertices that may be covered by new branch
  steminternal <- stem |> subset(!(x >= min(internal$x) & x <= max(internal$x) & y >= min(internal$y) & y <= max(internal$y))) |> subset(select=original)
  #remove tip of stem if branch too close to top
  if(max(subset(stem, type %in% c('base','mid','tip'))$y) - bht < 0.1){
    xdi <- mean(subset(steminternal, type %in% 'tip')$i)
    steminternal <- steminternal |> subset(!type %in% 'tip')
  }else{
    ylower2 <- max(subset(xd, y < min(newbase$y))$y)
    yupper2 <- min(subset(xd, y > max(newbase$y))$y)
    xdi <- mean(subset(xd, y %in% c(yupper2, ylower2))$i)
  }
  #identify where to insert new numbering sequence to maintain correct vertex order
  
  #assemble branch with new base, omitting internal vertices
  branchinternal <- branch |> subset(inside >= 0) |> subset(select=original) |> rbind(newbase[,original]) |> mutate(i = xdi + i/10000, type=paste0('b',type), inside = NULL, h=NULL,a=NULL)  |> arrange(i) 
  #append branch to stem with correct vertex order
  stemnew <- rbind(branchinternal,steminternal) |> arrange(i) 
  #renumber vertices
  stemnew <- mutate(stemnew, i=(1:nrow(stemnew)))
  return(stemnew)}

ggplot()+
  geom_polygon(data=steminternal, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=steminternal, aes(x=x, y=y), color='red')+
  geom_polygon(data=branch, aes(x=x, y=y), color='blue',fill='#00009950')+
  geom_point(data=branch, aes(x=x, y=y), color='blue')+
  geom_polygon(data=stemnew, aes(x=x, y=y), color='green',fill='#00990050')+
  geom_point(data=stemnew, aes(x=x, y=y), color='green')+
  coord_fixed()



#fractal tree with bend --- 
lth=4
wth = 0.5
sc = 0.6
stem <-  makeStem(lth,wth,0.2*wth,20)
stem <-  stem |> mutate(x = x + -0.2*cos(y/lth*2*pi), center = center + -0.2*cos(y/lth*2*pi))
branch <-  stem |> mutate(x=x*sc,y=y*sc,center=center*sc)
branch2 <-  stem |> mutate(x=x*sc*0.7,y=y*sc*0.7,center=center*sc*0.7)
branch3 <-  stem |> mutate(x=x*sc*0.3,y=y*sc*0.3,center=center*sc*0.3)

for(i in 1:3){
  tree <- attachBranch(stem, branch, -50, 2)
  branch <- branch2
  tree <- attachBranch(tree, branch, 50, 3)
  branch <- branch3
  #stem=tree
  tree <- attachBranch(tree, branch, -30, 4)
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




inc=8
s = (0:inc)/inc
i = 1:(inc+1)
angles <- s*2*pi/4

x = cos(angles)
y = sin(angles)
circle <- data.frame(i=i,x=cos(angles),y=sin(angles))
circle2 <- data.frame(i=i[(inc/2):inc],x=cos(angles[(inc/2):inc]),y=sin(angles[(inc/2):inc]))
triangle <- data.frame(i=i,x=-1*s+1,y=s*2)
verticle <- data.frame(i=i,x=1,y=s*2)
combo <- data.frame(i=i,x=(circle$x+triangle$x)/2, y=(circle$y+triangle$y)/2)
combo2 <- data.frame(i=i,x=(circle$x+verticle$x)/2, y=(circle$y+verticle$y)/2)
ggplot()+
  # geom_polygon(data=circle, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=circle, aes(x=x, y=y), color='red')+
  geom_point(data=triangle, aes(x=x, y=y), color='blue')+
  geom_point(data=verticle, aes(x=x, y=y), color='green')+
  geom_point(data=combo, aes(x=x, y=y), color='purple')+
  geom_point(data=circle2, aes(x=x, y=y), color='gold')+
  coord_fixed()

dfbranches <- data.frame(
  baseht = c(1.5,2,2.5,3,3.5,4),
  crownht = (circle[c(2:7),]$y)*2+3,
  crownwd = (circle[c(2:7),]$x))

dfbranches <- dfbranches |> mutate(branchlen = ((crownht - baseht)^2+crownwd^2)^0.5,
                                   angle=atan((crownht - baseht)/crownwd)*180/pi)


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



attachBranchsf <- function(stem, branch, angle, bht){
  #establish permanent columns to end up with
  original <- colnames(stem)
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
  #determine which stem vertices straddle the branch
  xd <- subset(stem, side %in% xside & type %in% c('tip','mid','base','bbase'))
  ylower <- max(subset(xd, y < min(branch[branch$type %in% 'base',]$y))$y)
  yupper <- min(subset(xd, y > max(branch[branch$type %in% 'base',]$y))$y)
  xupper <- subset(xd, y == yupper)$x
  xlower <- subset(xd, y == ylower)$x
  #shift branch to conform with twisted stem center position
  cshift <- mean(subset(xd, y == yupper | y == ylower)$center)
  branch <- branch |> mutate(x = x+cshift)
  
  #? bx <- (bht-ylower)/(yupper-ylower)*(xupper-xlower)+xlower
  #identify which branch vertices are inside stem
  branch <- branch |> mutate(inside = (x-((y-ylower)/(yupper-ylower)*(xupper-xlower)+xlower))*xsine)
  #identify vertices which straddle inside and outside of stem to find points of intersection
  internal <- branch |> mutate(near = inside^2, isinside = ifelse(inside < 0 | type %in% 'base', 'no','yes')) |> group_by(isinside, side) |> mutate(minnear = min(near)) |> ungroup() |> subset(near == minnear)
  #approximate location of branch stem intersection to establish new branch base
  newbase <- internal |> group_by(side) |> mutate(amt = 1/((inside - 0)/(max(inside)-min(inside)))^2) |> 
    summarise(x= sum(amt*x)/sum(amt), y= sum(amt*y)/sum(amt), i= mean(i), type='base', center=0, side=xside) 
  #remove stem vertices that may be covered by new branch
  steminternal <- stem |> subset(!(x >= min(internal$x) & x <= max(internal$x) & y >= min(internal$y) & y <= max(internal$y))) |> subset(select=original)
  #identify where to insert new numbering sequence to maintain correct vertex order
  ylower2 <- max(subset(xd, y < min(newbase$y))$y)
  yupper2 <- min(subset(xd, y > max(newbase$y))$y)
  xdi <- mean(subset(xd, y %in% c(yupper2, ylower2))$i)
  #assemble branch with new base, omitting internal vertices
  branchinternal <- branch |> subset(inside >= 0) |> subset(select=original) |> rbind(newbase[,original]) |> mutate(i = xdi + i/10000, type=paste0('b',type), inside = NULL, h=NULL,a=NULL)  |> arrange(i) 
  #append branch to stem with correct vertex order
  stemnew <- rbind(branchinternal,steminternal) |> arrange(i) 
  #renumber vertices
  stemnew <- mutate(stemnew, i=(1:nrow(stemnew)))
  return(stemnew)}










tree2 <- data.frame(x=st_coordinates(polygon)[,1],y=st_coordinates(polygon)[,2])
ggplot()+
  geom_polygon(data=tree2, aes(x=x, y=y), color='green',fill='#00990050')+
  coord_fixed()