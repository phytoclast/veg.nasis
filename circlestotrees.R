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
                   center = 0)
  df[c(1,nrow(df)),]$type <- 'base'
  df[c(nrow(df)/2,nrow(df)/2+1),]$type <- 'tip'
  return(df)
}
angle=-50
bht=2.5

attachBranch <- function(stem, branch, angle, bht){
angle = angle/360*2*pi
xside <- ifelse(angle>0,'R','L')
xsine <- ifelse(angle>0,1,-1)
branch <- branch |> mutate(h=(x^2+y^2)^0.5,a = acos(y/h),a=ifelse(x < 0,-1*a,a),
                      x = h*sin(a+angle), y = h*cos(a+angle))

branch <- branch |> mutate(y = y+bht)

xd <- subset(stem, side %in% xside)
ylower <- max(subset(xd, y < min(branch[branch$type %in% 'base',]$y))$y)
yupper <- min(subset(xd, y > max(branch[branch$type %in% 'base',]$y))$y)
xupper <- subset(xd, y == yupper)$x
xlower <- subset(xd, y == ylower)$x
cshift <- mean(subset(xd, y == yupper | y == ylower)$center)

bx <- (bht-ylower)/(yupper-ylower)*(xupper-xlower)+xlower
branch <- branch |> mutate(x = x+cshift)

branch <- branch |> mutate(inside = (x-((y-ylower)/(yupper-ylower)*(xupper-xlower)+xlower))*xsine)
#remove internal vertices

internal <- branch |> mutate(near = inside^2, isinside = ifelse(inside > 0, 'no','yes')) |> group_by(isinside, side) |> mutate(minnear = min(near)) |> ungroup() |> subset(near == minnear)

newbase <- internal |> group_by(side) |> mutate(amt = 1/((inside - 0)/(max(inside)-min(inside)))^2) |> 
  summarise(x= sum(amt*x)/sum(amt), y= sum(amt*y)/sum(amt), i= mean(i), type='base', center=0) 

steminternal <- stem |> subset(!(x >= min(internal$x) & x <= max(internal$x) & y >= min(internal$y) & y <= max(internal$y)))
xdi <- mean(subset(xd, y %in% c(yupper, ylower))$i)

branchinternal <- branch |> subset(inside >= 0) |> subset(select=colnames(stem)) |> rbind(newbase[,colnames(stem)]) |> mutate(i = xdi + i/1000, type=paste0('b',type), inside = NULL, h=NULL,a=NULL) 

stemnew <- rbind(branchinternal,steminternal) |> arrange(i) 
stemnew <- mutate(stemnew, i=(1:nrow(stemnew)))
return(stemnew)}


stem <-  makeStem(4,1,0.1,10)
branch <- makeStem(2,.2,0.01,10)
branch2 <- makeStem(1.5,.2,0.01,10)
tree <- attachBranch(stem, branch, -90, 2)
tree <- attachBranch(tree, branch2, 85, 3)
tree2 <- tree |> mutate(x=x/2, y=y/2)
tree <- attachBranch(stem, tree2, -100, 2)
tree <- attachBranch(tree, tree2, 90, 3)

ggplot()+
  geom_polygon(data=stem, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=branch, aes(x=x, y=y), color='blue',fill='#00009950')+
  geom_point(data=branch, aes(x=x, y=y), color='blue')+
  geom_polygon(data=tree, aes(x=x, y=y), color='green',fill='#00990050')+
  geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()

lth=4
wth = 0.5
sc = .25
stem <-  makeStem(lth,wth,0.1*wth,20)
stem <-  stem |> mutate(x = x + -0.2*cos(y/lth*2*pi), center = center + -0.2*cos(y/lth*2*pi))
branch <-  stem |> mutate(x=x*sc,y=y*sc,center=center*sc)
branch2 <-  stem |> mutate(x=x*sc*0.7,y=y*sc*0.7,center=center*sc*0.7)

for(i in 1:4){
  tree <- attachBranch(stem, branch, -50, 2)
  tree <- attachBranch(tree, branch2, 50, 3)
  branch <- tree |> mutate(x=x*sc,y=y*sc,center=center*sc)
  branch2 <- tree |> mutate(x=x*sc*.7,y=y*sc*.7,center=center*sc*.7)
}

ggplot()+
  # geom_polygon(data=stem, aes(x=x, y=y), color='red',fill='#99000050')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=tree, aes(x=x, y=y), color='green',fill='#00990050')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()