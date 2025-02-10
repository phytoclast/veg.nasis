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


#add options to add branch constraint by height and crown width rather than just angle

attachBranch <- function(stem, branch, angle=90, bht, tht=NA, tx=NA){
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
  # if(bwd > swd) {branch <- branch |> mutate(x = x*swd/bwd)}
  
  #convert angle to radians
  angle = angle/360*2*pi  
  
  #optional establish branch angle based on branch distance from stem
  if(!is.na(tht) & !is.na(tx)){
  l = ((tht-bht)^2+(tx)^2)^0.5
  bmaxlen <- max((branch$x^2+branch$y^2)^0.5)
  #resize branch to fit specified height and width
  branch <- branch |> mutate(x=x*l/bmaxlen, y=y*l/bmaxlen)
  angle <- asin(tx/l)
  }
  

  
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


skewStem <- function(stem, amp=0.2, phase=0, waves=1){
  maxstem <- max(stem$y)
  minstem <- min(stem$y)
  lth <- maxstem - minstem
  stem <-  stem |> mutate(x = x + amp*cos((y/lth+phase)*2*pi*waves), 
                          center = center + amp*cos((y/lth+phase)*2*pi*waves))
  return(stem)
}

makeCrownShape <- function(ht.max=5, ht.min=1, crwd=2, dbh, crshape=c('pyramid','dome','round','column'), n=5, bu=0.8, bl=0, opposite = FALSE){
  
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
  geom_polygon(data=tree, aes(x=x, y=y), color='brown',fill='#99500090')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=crown, aes(x=x, y=y), color='green',fill='#00990090')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()

#------conifer

crshape = c('pyramid','dome','round','column')

crshape = c('pyramid','column')
shapes <- makeCrownShape(ht.max=10, ht.min=3, crwd=3, dbh=0.5, n=7, bu=1, bl=0.2, crshape=crshape,
                        opposite = T)
shapes <- subset(shapes, !(a > 175 | a < -175 | a == 0)  & l> 0.3)
stem <-  makeStem(10,0.5,0.01,10)
for(i in 1:nrow(shapes)){#i=1
  branch <- makeStem(shapes$l[i], shapes$d[i]*.5,0.01,10)
  bpos <- max(branch$y)
  branch2 <- branch |> mutate(x=x*0.3,y=y*0.3)
  branch3 <- branch |> mutate(x=x*0.2,y=y*0.2)
  
  branch <- skewStem(branch, amp=ifelse(shapes$a[i] >= 0,-0.07*(shapes$s[i]*-1+1),0.07*(shapes$s[i]*-1+1)), 
                     phase=0, waves=1)
  branchA <- attachBranch(branch, branch2, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.3)
  branchA <- attachBranch(branchA, branch3, ifelse(shapes$a[i] >= 0,-30,30), bpos*0.7)
  
  
  stem <- attachBranch(stem, branchA, shapes$a[i], shapes$by[i])
}
crown <- stem |> subset(grepl('tip',type) | type %in% 'bbase')


ggplot()+
  geom_polygon(data=stem, aes(x=x, y=y), color='brown',fill='#99500090')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=crown, aes(x=x, y=y), color='green',fill='#00990090')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()

#------cottonwood

crshape = c('pyramid','dome','round','column')

crshape = c('dome')
shapes <- makeCrownShape(ht.max=10, ht.min=5, crwd=6, dbh=0.5, n=3, bu=0.5, bl=0, crshape=crshape,
                        opposite = F)
shapes <- subset(shapes, !(a > 175 | a < -175 | a==0)  & l> 0.3)
stem <-  makeStem(10,0.5,0.01,25)
stem <- skewStem(stem, amp=-0, 
                   phase=0, waves=0.5)

for(i in 1:nrow(shapes)){#i=1
  branch <- makeStem(shapes$l[i], shapes$d[i]*.5,0.01,15)
  bpos <- max(branch$y)
  branch2 <- branch |> mutate(x=x*0.3,y=y*0.3)
  branch3 <- branch |> mutate(x=x*0.2,y=y*0.2)
  
  branch <- skewStem(branch, amp=ifelse(shapes$a[i] >= 0,-0.2,0.2), 
                     phase=0, waves=1)
  branchA <- attachBranch(branch, branch2, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.5)
  branchA <- attachBranch(branchA, branch3, ifelse(shapes$a[i] >= 0,-30,30), bpos*0.6)
  
  
  stem <- attachBranch(stem, branchA, shapes$a[i], shapes$by[i], tht = shapes$ty[i], tx = shapes$tx[i])
}

branchB <- stem |> mutate(x=x*1,y=y) |> skewStem(amp=0.3)

stem2 <-  makeStem(5,0.5,0.3,30)
tree <- attachBranch(stem2, branchB, bht=4.5, tht = 11, tx=-4)
tree <- attachBranch(tree, branchB, bht=5, tht = 11, tx=3)

crown <- tree |> subset(grepl('tip',type))
crown2 <- crown[chull(x=crown$x, y=crown$y),]

ggplot()+
  geom_polygon(data=tree, aes(x=x, y=y), color='brown',fill='#99500050')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=crown2, aes(x=x, y=y), color='green',fill='#00990050')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()






x = c(-2,-2,0,2,2)
y = c(0,1,0.5,1,0)
df <- data.frame(x=x,y=y)

#Convex hull (trying to modify to allow some limited concavity)
#eliminate redundant points via rounding
df <- df |> mutate(x = round(x,2),y = round(y,2)) 
x = df$x 
y = df$y 

#set up data frame for calculations 
df <-  data.frame(a=NA,l1=NA,l2=NA,x=x,y=y) |> unique()
#track original order for troubleshooting
df <- df |> mutate(q = 1:nrow(df))
rm(x,y)

#find starting point at bottom of plot
miny <- min(df$y)
minx <- min(df[df$y==miny,]$x)
df <- mutate(df, s = ifelse(x==minx & y==miny, 1,NA))
i=1
refx1 <- df[df$s %in% i,]$x
refy1 <- df[df$s %in% i,]$y

df <- df |> mutate(l1 = ((x-refx1)^2+(y-refy1)^2)^0.5,
                   a = asin((x-refx1)/l1), a=ifelse(l1==0,0,a),
                   deg = a/2/pi*360)
amax <- max(df[!df$s %in% c(i,i-1),]$a)
amin <- min(df[!df$s %in% c(i,i-1),]$a)
nxtpt <- df |> subset((a == amax | a == amin) & is.na(s))
minl <- min(nxtpt$l1)
nxtpt <- nxtpt |> subset(l1 == minl & is.na(s))

df <- df |> mutate(s = ifelse(x == nxtpt$x & y == nxtpt$y, i+1, s))
#establish variable that will indicate when to stop
check <- TRUE 

#loop through subsequent points to identify which points to retain
for (i in 2:nrow(df)){#   i=4
  if(check){
    #current and previous points
    refx1 <- df[df$s %in% (i-1),]$x
    refy1 <- df[df$s %in% (i-1),]$y
    refx2 <- df[df$s %in% i,]$x
    refy2 <- df[df$s %in% i,]$y
    
    #distance between current and previous points
    l3 = ((refx1-refx2)^2+(refy1-refy2)^2)^0.5
    #angle of first seqment
    a0=acos((refx1-refx2)/l3)
    a0=ifelse(refy1-refy2 >=0,a0,-1*a0)
    #distance between current and previous points to every other point to get the total angle of current point relative to potential next point
    df <- df |> mutate(l1 = ((x-refx1)^2+(y-refy1)^2)^0.5,
                       l2 = ((x-refx2)^2+(y-refy2)^2)^0.5,
                       a=acos((l2^2+l3^2-l1^2)/(2*l2*l3)),
                       a=ifelse(is.na(a),0,a),
                       anew=cos((x-refx2)/l2),
                       anew=ifelse(y-refy2 >=0,anew,-1*anew),
                       adif=anew-(a0-pi),#angle difference between 2nd and 1st segments to see if left or right
                       deg = a/2/pi*360,
                       a1 = ifelse(l2 > 3,a/2,a),
                       a1 = ifelse(l2 > 3,adif/2,adif))#experimental -  trying to find a alternative path based on distance 
    #identify which angle is the largest to ensure convexity
    amax <- max(subset(df, !s %in% c(i,i-1))$a1)
    #set trigger if next point is already assigned, which means perimeter has been closed
    check <- is.na(subset(df,a1 == amax)$s)
    #assign next point number
    df <- df |> mutate(s = ifelse(a1 == amax & is.na(s), i+1, s))
  }}
#exclude excess points and sort the by path order
df2 <- subset(df,!is.na(s)) |> arrange(s) 
ggplot()+
  geom_polygon(data=df, aes(x=x, y=y), color='blue',fill='#50500050')+
  geom_point(data=df, aes(x=x, y=y), color='blue',fill='#50500050')+
  geom_polygon(data=df2, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=df2, aes(x=x, y=y), color='red',fill='#99000050')+
  coord_fixed()






############XXXXXXXXXXXXXXXXXXXXXxx trial 2


x = c(-1,-1.5,0,1,3,1,0.5)
y = c(0,1,0.5,1,0.5,0,0.2)
df <- data.frame(x=x,y=y)

df <- df |> mutate(q = 1:nrow(df))
rm(x,y)
check <- TRUE 
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
for(i in 1:2){#nrow(df) i=1
  if(check){
    x0 <- df[df$s %in% (i-1),]$x
    y0 <- df[df$s %in% (i-1),]$y
    x1 <- df[df$s %in% i,]$x
    y1 <- df[df$s %in% i,]$y
    l0 = ((x1-x0)^2+(y1-y0)^2)^0.5
    a0 = acos((x1-x0)/l0)/2/pi*360
    a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
    df <- df |> mutate(l1 = ((x-x1)^2+(y-y1)^2)^0.5,
                       a1=acos((x-x1)/l1)/2/pi*360,
                       a1=ifelse(y-y1 >=0,a1,-1*a1),
                       a1= a1-a0,
                       a1= ifelse(a1 > 180, 360-a1,ifelse(a1 < -180, -360-a1, a1)),
                       a1= case_when(x >= x1 & y >= y1 & y1 < y0 ~ -a1,
                                     x < x1 & y < y1 & y1 >= y0 ~ -a1,
                                     TRUE ~ a1))
    
    
    amax = max(subset(df, !s %in% c(i-1,i) )$a1)#& a1 < 90
    lmin = min(subset(df, !s %in% c(i-1,i) & a1 %in% amax)$l1)
    df <- df |> mutate(s = ifelse(s %in% 0,NA,s))
    check <- is.na(subset(df,a1 %in% amax & l1 %in% lmin)$s)
    df <- df |> mutate(s = ifelse(a1 == amax & is.na(s) & l1 %in% lmin, i+1, s))
  }
}

df2 <- subset(df,!is.na(s)) |> arrange(s) 
ggplot()+
  geom_polygon(data=df, aes(x=x, y=y), color='blue',fill='#20200060')+
  geom_point(data=df, aes(x=x, y=y), color='blue',fill='#20200060')+
  geom_polygon(data=df2, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=df2, aes(x=x, y=y), color='red',fill='#99000050')+
  coord_fixed()







#trial 3 -----
x = c(-1,-1.5,0,1,3,1,0.5)
y = c(0,1,0.51,1,0.5,0,0.2)
n=3
df <- data.frame(x=tree$x,y=tree$y)
#df <- data.frame(x=x,y=y)

df <- df |> mutate(q = 1:nrow(df))
rm(x,y)

#convex hull ----
check <- TRUE 
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
for(i in 1:nrow(df)){#nrow(df) i=2
  if(check){
    x0 <- df[df$s %in% (i-1),]$x
    y0 <- df[df$s %in% (i-1),]$y
    x1 <- df[df$s %in% i,]$x
    y1 <- df[df$s %in% i,]$y
    l0 = ((x1-x0)^2+(y1-y0)^2)^0.5
    a0 = acos((x1-x0)/l0)
    a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
    df <- df |> mutate(
      l2 = ((x-x0)^2+(y-y0)^2)^0.5,
      l1 = ((x-x1)^2+(y-y1)^2)^0.5,
      a1=acos((l1^2+l0^2-l2^2)/(2*l1*l0)),
      a1=ifelse(is.na(a1),0,a1))
    
    
    amax = max(subset(df, !s %in% c(i-1,i) )$a1)
    lmin = min(subset(df, !s %in% c(i-1,i) & a1 %in% amax)$l1)
    df <- df |> mutate(s = ifelse(s %in% 0,NA,s))
    check <- is.na(subset(df,a1 %in% amax & l1 %in% lmin)$s)
    df <- df |> mutate(s = ifelse(a1 == amax & is.na(s) & l1 %in% lmin, i+1, s))
  }
}

#concave hull first degree ----
df <- df |> mutate(s1 = s)

smax <- max(df$s, na.rm = TRUE)

for(i in 1:smax){
  i0=ifelse(i == 1,smax,i-1)
  x0 <- df[df$s %in% (i0),]$x
  y0 <- df[df$s %in% (i0),]$y
  x1 <- df[df$s %in% i,]$x
  y1 <- df[df$s %in% i,]$y
  l0 <- ((x1-x0)^2+(y1-y0)^2)^0.5
  a0 = acos((x1-x0)/l0)
  a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
  df <- df |> mutate(xr= x-x0,
                     yr= y-y0,
                     h=(xr^2+yr^2)^0.5,
                     a1=acos(yr/h),
                     a1=ifelse(xr >=0,a1,-1*a1),
                     a1= a1+a0,
                     xr = ifelse(h==0,0,h*sin(a1)), 
                     yr = ifelse(h==0,0,h*cos(a1)),
                     xr = xr/l0)
  bottom <- min(l0,abs(min(df[df$xr > -1.25 & df$xr < 1.25 ,]$yr)))*-1
  df <- df |> mutate(yr= yr/(bottom))
  
  for(j in 1:n){
    psmin1 <- (subset(df, xr >= (j-1)/n & xr < j/n)$yr)
    
    smin1 <- ifelse(length(psmin1)>0,min(psmin1),1)
    
    df <- df |> mutate(s1 = case_when(xr >= (j-1)/n & xr < j/n & yr %in% smin1 & yr < 0.5 ~ i0+j/n/100,
                                      TRUE ~ s1))
  }
}
#concave hull second degree ----
n=5
df$s2 <- df$s1
gather <- df$s1
notblank <- !is.na(gather)
justthese <- gather[notblank]
smax <- length(justthese)
for(k in 1:smax){
  i <- justthese[k]
  i0=ifelse(k == 1,justthese[length(justthese)],justthese[k-1])
  x0 <- df[df$s1 %in% (i0),]$x
  y0 <- df[df$s1 %in% (i0),]$y
  x1 <- df[df$s1 %in% i,]$x
  y1 <- df[df$s1 %in% i,]$y
  l0 <- ((x1-x0)^2+(y1-y0)^2)^0.5
  a0 = acos((x1-x0)/l0)
  a0 = ifelse(y1 - y0 >=0,a0,-1*a0)
  df <- df |> mutate(xr= x-x0,
                     yr= y-y0,
                     h=(xr^2+yr^2)^0.5,
                     a1=acos(yr/h),
                     a1=ifelse(xr >=0,a1,-1*a1),
                     a1= a1+a0,
                     xr = ifelse(h==0,0,h*sin(a1)), 
                     yr = ifelse(h==0,0,h*cos(a1)),
                     xr = xr/l0)
  bottom <- min(l0,abs(min(df[df$xr > -1.25 & df$xr < 1.25 ,]$yr)))*-1
  df <- df |> mutate(yr= yr/(bottom))
  j=1
  for(j in 1:n){
    psmin1 <- (subset(df, xr >= (j-1)/n & xr < j/n)$yr)
    
    smin1 <- ifelse(length(psmin1)>0,min(psmin1),1)
    
    df <- df |> mutate(s2 = case_when(xr >= (j-1)/n & xr < j/n & yr %in% smin1 & yr < 0.5 ~ i0+j/n/10000,
                                      TRUE ~ s2))
  }
}

x <- df$s1
xna <- !is.na(x)
xnew <- order(x)
xnao <- xna[xnew]
ret <- xnew[xnao]



df2 <- subset(df,!is.na(s1)) |> arrange(s1) 
df3 <- df[ret,]
ggplot()+
  geom_polygon(data=df, aes(x=x, y=y), color='blue',fill='#20200060')+
  geom_point(data=df, aes(x=x, y=y), color='blue',fill='#20200060')+
  geom_polygon(data=df2, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_point(data=df2, aes(x=x, y=y), color='red',fill='#99000050')+
  geom_polygon(data=df3, aes(x=x, y=y), color='green',fill='#00990050')+
  geom_point(data=df3, aes(x=x, y=y), color='green',fill='#00990050')+
  coord_fixed()

ggplot()+
  geom_polygon(data=df, aes(x=x, y=y), color='blue',fill='#20200060')+
  geom_point(data=df, aes(x=x, y=y), color='blue',fill='#20200060')+
  geom_polygon(data=df, aes(x=xr, y=yr), color='red',fill='#99000050')+
  geom_point(data=df, aes(x=xr, y=yr), color='red',fill='#99000050')+
  coord_fixed()

x <- c(1,2,9,NA,3,NA,8,7,6)
xna <- !is.na(x)
xnew <- order(x)
xnao <- xna[xnew]
return(xnew[xnao])
xfinal <-  x[xnew[xnao]]
xfinal