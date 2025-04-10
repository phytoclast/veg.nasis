library(vegnasis)
library(ggplot2)
shps <-  vegnasis::shapes

tree1 <- make_tree(ht.max=10, ht.min=5, crwd=5, dbh=20, crshape='blob1', stshape='trunk')



ggplot()+
  geom_polygon(data=subset(tree1, obj %in%  'stem'), aes(x=x, y=z, color=color,fill=fill))+
  geom_polygon(data=subset(tree1, obj %in%  'crown'), aes(x=x, y=z, color=color,fill=fill))+
  scale_color_manual(values=unique(tree1$color),breaks=unique(tree1$color))+
  scale_fill_manual(values=unique(tree1$fill),breaks=unique(tree1$fill))+
  theme(legend.position = "none")+
  coord_fixed()

##############################
crshape = c('pyramid','dome','round','column')
ht.max=15
ht.min=5
crwd=5
dbh=40
bu=1
bl=0.1
opposite=T
oppfactor = ifelse(opposite,1,2)
n = pmax(floor(10*(ht.max-ht.min)*(bu-bl)/10*oppfactor/crwd*5),1)
crshape = c('pyramid')

bf <- ifelse(opposite, 5/n,10/n)
shapes <- makeCrownShape(ht.max=ht.max,ht.min=ht.min, crwd=crwd, dbh=dbh/100, n=n, bu=bu, bl=bl, crshape=crshape,opposite = opposite)

shapes <- subset(shapes, !(a > 175 | a < -175 | a == 0)  & l> 0.15 & by < ht.max)
stem <-  makeStem(ht.max,dbh/100,0.05,15)
cstem <- stem
for(i in 1:nrow(shapes)){#i=1
  branch <- makeStem(shapes$l[i], shapes$d[i]*0.5,0.01,10)
  bpos <- max(branch$y)
  twigA <- branch |> mutate(x=x*0.4*bf,y=y*0.4*bf)
  twigB <- branch |> mutate(x=x*0.3*bf,y=y*0.3*bf)
  twigC <- branch |> mutate(x=x*0.2*bf,y=y*0.2*bf)
  branchx1 <- branch |> mutate(x=x*0.3*bf,y=y*0.3*bf)
  branch <- skewStem(branch, amp=ifelse(shapes$a[i] >= 0,-0.07*(shapes$s[i]*-1+1),0.07*(shapes$s[i]*-1+1)),
                     phase=0, waves=1)
  # branchA <- attachBranch(branch,  twigB, ifelse(shapes$a[i] >= 0,40,-40), bpos*0.15)#lower
  # branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[i] >= 0,40,-40), bpos*0.3)#lower
  branchA <- attachBranch(branch, twigB, ifelse(shapes$a[i] >= 0,35,-35), bpos*0.4)#lower
  # branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[i] >= 0,35,-35), bpos*0.5)#lower
  # branchA <- attachBranch(branchA, twigB, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.6)#lower
  # branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.7)#lower
  branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[i] >= 0,-30,30), bpos*0.7)#upper
  # branchA <- attachBranch(branchA, twigB, ifelse(shapes$a[i] >= 0,-30,30), bpos*0.5)#upper
  #basal branch twigs to attach crown
  branchB <- attachBranch(branchA, branchx1, ifelse(shapes$a[i] >= 0,-90,90), bpos*0.05)
  branchB <- attachBranch(branchB, branchx1, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.05)

  stem <- attachBranch(stem, branchA, shapes$a[i], shapes$by[i])#branches to show bare
  cstem <- attachBranch(cstem, branchB, shapes$a[i], shapes$by[i])#branches to attach crown
}
crown <- cstem |> subset(grepl('tip',type))
crown2 <- cavhull(x=crown$x,y=crown$y, concavity = 1, curvy = T, mag = 1, minspan = 0.0, maxdepth = 1)
saveRDS(crown,'C:/scripts/veg.nasis/df.rds')
ggplot()+
  geom_polygon(data=stem, aes(x=x, y=y), color='brown',fill='#99500090')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  # geom_polygon(data=crown2, aes(x=x, y=y), color='green',fill='#00990090')+
  geom_polygon(data=crown2, aes(x=x, y=y), color='green',fill='#00990090')+
  geom_point(data=crown, aes(x=x, y=y), color='green')+
  coord_fixed()



###########################################
crshape = c('pyramid','dome','round','column')
ht.max=15
ht.min=5
crwd=10
dbh=25
bu=0.5#as crown widens, make top and bottom lower
bl=-0.3#make sure tree is tall enough to support negative
opposite=F
oppfactor = ifelse(opposite,1,2)
n = pmax(floor(7*(ht.max-ht.min)*(bu-bl)/10*oppfactor/crwd*5),1)
crshape = c('dome')
ca <- pmin(1,(crwd*(ht.max-ht.min)*(bu-bl)/80))

shapes <- makeCrownShape(ht.max=ht.max,ht.min=ht.min, crwd=crwd, dbh=dbh/100, n=n, bu=bu, bl=bl, crshape=crshape,opposite = opposite)
shapes <- subset(shapes, !(a > 175 | a < -175 | a == 0)  & l> 0.15 & by < ht.max)

bf <- ifelse(opposite, 5/n,10/n)#twig size based on number of branches
af <- pmin(1,mean(abs(shapes$a))/30)#twig size based on angles

stem <-  makeStem(ht.max,dbh/100,0.02,15)
stem <- skewStem(stem, amp = 0.1, phase = 0, waves = 1)

branch <- makeStem(crwd*0.1*ca*af, dbh/100*0.1,0.01,10)#small branch near top
stem <- attachBranch(stem, branch, -30*af, ht.max*0.93)
cstem <- stem
for(i in 1:nrow(shapes)){#i=1
  branch <- makeStem(shapes$l[i], shapes$d[i]*0.5,0.01,10)
  bpos <- max(branch$y)
  twigA <- branch |> mutate(x=x*0.4*bf,y=y*0.4*bf*af)
  twigB <- branch |> mutate(x=x*0.2*bf,y=y*0.2*bf*af)
  twigB <- skewStem(twigB, amp=ifelse(shapes$a[i] >= 0,
                                      -0.05*(af*bf*-1+1),
                                      0.05*(af*bf*-1+1)),
                    phase=0.5, waves=1.5)
  twigC <- branch |> mutate(x=x*0.1*bf,y=y*0.1*bf*af)
  twigC <- skewStem(twigC, amp=ifelse(shapes$a[i] >= 0,
                                      -0.03*(af*bf*-1+1),
                                      0.03*(af*bf*-1+1)),
                    phase=0.5, waves=1.5)
  branchx1 <- branch |> mutate(x=x*0.3*bf,y=y*0.3*bf*af)
  branch <- skewStem(branch, amp=ifelse(shapes$a[i] >= 0,-0.07*(shapes$s[i]*-1+1),0.07*(shapes$s[i]*-1+1)), phase=0.5, waves=2)
  branchA <- attachBranch(branch, twigB, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.3)#lower
  branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[i] >= 0,-30,30), bpos*0.7)#upper
  #basal branch twigs to attach crown
  branchB <- attachBranch(branchA, branchx1, ifelse(shapes$a[i] >= 0,-90,90), bpos*0.05)
  branchB <- attachBranch(branchB, branchx1, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.05)
  stem <- attachBranch(stem, branchA, shapes$a[i], shapes$by[i])#branches to show bare
  cstem <- attachBranch(cstem, branchB, shapes$a[i], shapes$by[i])#branches to attach crown
}
crown <- cstem |> subset(grepl('tip',type))
circles <- data.frame(a=(0:7)/7*2*pi) |> mutate(cx=cos(a),cy=sin(a))
crown2 <- merge(crown, circles)# |> mutate(x=x+0.2*cx,y=y+0.2*cy)
crown2 <- cavhull(x=crown2$x,y=crown2$y, concavity = 1, curvy = T, maxdepth = 1, minspan = 2, 
                   mag = 3)

ggplot()+
  geom_polygon(data=stem, aes(x=x, y=y), color='brown',fill='#99500090')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
  geom_polygon(data=crown2, aes(x=x, y=y), color='green',fill='#00990090')+
  # geom_polygon(data=crown, aes(x=x, y=y), color='green',fill='#00990090')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()



veg.raw <-  vegnasis::nasis.veg
veg <- clean.veg(veg.raw)

veg.select <- subset(veg,  grepl('2022MI165023.P',plot))
plants <- grow_plants(veg.select)

veg_profile_plot(plants, unit='m',  skycolor = 'white', fadecolor = 'lightgray', gridalpha = 0.1, groundcolor = 'darkgray')


veg.select <- subset(veg,  grepl('2022MI165021.P',plot))

taxon <- c('Acer rubrum', 'Pinus resinosa')
crfill <- c("#80991A01","#80991A")
stfill <- c('gray',"#B36666")
crshape <- c('hardwood','boreal')
override <- data.frame(taxon=taxon,stfill=stfill,crfill=crfill,crshape=crshape)
veg.select <- veg.select |> left_join(override)

plants <- grow_plants(veg.select)
plants2 <- plants |> subset(!(shape %in% 'hardwoodcrown' & obj %in% 'crown'))
veg_profile_plot1(plants2)

###############################
library(ggplot2)
library(vegnasis)
newtree <-  tree.001a(ht.max=15, ht.min=3, crwd=5, dbh=35)

ggplot()+
  geom_polygon(data=subset(newtree, obj %in% 'stem'), aes(x=x, y=z), color='brown',fill='#99500090')+
  geom_polygon(data=subset(newtree, obj %in% 'crown'), aes(x=x, y=z), color='green',fill='#00990090')+
  coord_fixed()


###################################
tree.001a <- function(ht.max,
                      ht.min,
                      crwd,
                      dbh){
  ht.max=15
  ht.min=5
  crwd=5
  dbh=35
  tip=0.05
  bu=1
  bl=0.1
  opposite=T
  oppfactor = ifelse(opposite,1,2)
  n = pmax(floor(20*(ht.max-ht.min)*(bu-bl)/10*oppfactor/crwd*5),1)
  crshape = c('pyramid','dome')

  bf <- ifelse(opposite, 5/n,10/n)
  shapes <- makeCrownShape(ht.max=ht.max,ht.min=ht.min, crwd=crwd, dbh=dbh/100, tip=tip, n=n, bu=bu, bl=bl, crshape=crshape,opposite = opposite)

  #alternate branch length
  q<-c(1,1,0.5,0.5)
  nq <- nrow(shapes)
  qq <- rep(q,nq)[1:nq]
  shapes <- shapes |> mutate(l=qq*l)
  shapes <- shapes |> mutate(tx2=bx+cos((90-a)/360*2*pi)*l,
                             ty2=by+sin((90-a)/360*2*pi)*l)

  #filter bad branches
  shapes <- subset(shapes, !(a > 175 | a < -175 | a == 0)  & l>= 0.2 & by < ht.max-d*2)
  stem <-  makeStem(ht.max,dbh/100,tip,50)
  cstem <- stem
  for(j in 1:nrow(shapes)){#  nrow(shapes)    j=29
    branch <- makeStem(shapes$l[j], shapes$d[j]*0.5,0.01,10)
    bpos <- max(branch$y)
    twigA <- branch |> mutate(x=x*0.4*bf,y=y*0.4*bf)
    twigB <- branch |> mutate(x=x*0.4*bf,y=y*0.4*bf)
    twigC <- branch |> mutate(x=x*0.3*bf,y=y*0.3*bf)
    branchx1 <- branch |> mutate(x=x*0.3*bf,y=y*0.3*bf)
    branch <- skewStem(branch, amp=ifelse(shapes$a[j] >= 0,-0.07*(shapes$s[j]*-1+1),0.07*(shapes$s[j]*-1+1)), phase=0, waves=1)
    # branchA <- attachBranch(branch,  twigB, ifelse(shapes$a[j] >= 0,40,-40), bpos*0.15)#lower
    # branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[j] >= 0,40,-40), bpos*0.3)#lower
    branchA <- attachBranch(branch, twigB, ifelse(shapes$a[j] >= 0,35,-35), bpos*0.4)#lower
    # branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[j] >= 0,35,-35), bpos*0.5)#lower
    # branchA <- attachBranch(branchA, twigB, ifelse(shapes$a[j] >= 0,30,-30), bpos*0.6)#lower
    # branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[j] >= 0,30,-30), bpos*0.7)#lower
    branchA <- attachBranch(branchA, twigC, ifelse(shapes$a[j] >= 0,-30,30), bpos*0.7)#upper
    # branchA <- attachBranch(branchA, twigB, ifelse(shapes$a[j] >= 0,-30,30), bpos*0.5)#upper

    # branchA <- branch
    #basal branch twigs to attach crown
    branchB <- attachBranch(branchA, branchx1, ifelse(shapes$a[j] >= 0,-90,90), bpos*0.05)
    branchB <- attachBranch(branchB, branchx1, ifelse(shapes$a[j] >= 0,30,-30), bpos*0.05)
    branchA <- mutate(branchA, type = paste0(type,j))
    # branch <- mutate(branch, type = paste0(type,j))
    if(shapes$l[j] < .1){
      stem <- attachBranch(stem, branch, angle=shapes$a[j], bht=shapes$by[j])#branches to show bare
      cstem <- attachBranch(cstem, branch, shapes$a[j], shapes$by[j])#branches to attach crown
    }else{
      stem <- attachBranch(stem, branchA, shapes$a[j], shapes$by[j])#branches to show bare
      cstem <- attachBranch(cstem, branchB, shapes$a[j], shapes$by[j])#branches to attach crown
    }
  }
  crown <- cstem |> subset(grepl('tip',type))
  rownames(crown) <- crown$i
  rownames(stem) <- stem$i

  crown <- crown |> mutate(z=y, shape = 'borealcrown',  fill='green', color='darkgreen', obj='crown', ptord=i) |> select(c("x","z","shape","fill","color","obj","ptord"))
  stem <- stem |> mutate(z=y, shape = 'borealstem',  fill='orange', color='brown', obj='stem', ptord=i) |> select(c("x","z","shape","fill","color","obj","ptord"))
  newtree = rbind(stem,crown)
  return(newtree)
}
###############################
library(ggplot2)
library(vegnasis)
newtree <-  tree.001a(ht.max=15, ht.min=3, crwd=5, dbh=35)

ggplot()+
  geom_polygon(data=subset(newtree, obj %in% 'stem'), aes(x=x, y=z), color='brown',fill='#99500090')+
  geom_polygon(data=subset(newtree, obj %in% 'crown'), aes(x=x, y=z), color='green',fill='#00990090')+
  coord_fixed()

##########
angle=shapes$a[j]; bht=shapes$by[j]; tht = NA; tx = NA

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
  bbase <- subset(branch, grepl('base',type))
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


  #determine if stem if branch too close to top, less than branch width
  if(stemax - bht < bwd*2){
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



  ggplot()+
    geom_polygon(data=stemnew, aes(x=x, y=y), color='brown',fill='#99500090')+
    coord_fixed()

