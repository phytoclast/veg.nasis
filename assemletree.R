library(vegnasis)
library(ggplot2)
shps <-  vegnasis::shapes

tree1 <- make_tree(ht.max=10, ht.min=5, crwd=5, dbh=20, crshape='blob1', stshape='trunk')

vegnasis::makeCrownShape()


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
dbh=25
bu=1
bl=0.1
opposite=T
oppfactor = ifelse(opposite,1,2)
n = pmax(floor(10*(ht.max-ht.min)*(bu-bl)/10*oppfactor/crwd*5),1)
crshape = c('pyramid')

bf <- ifelse(opposite, 5/n,10/n)
shapes <- makeCrownShape(ht.max=ht.max,ht.min=ht.min, crwd=crwd, dbh=dbh/100, n=n, bu=bu, bl=bl, crshape=crshape,opposite = opposite)

shapes <- subset(shapes, !(a > 175 | a < -175 | a == 0)  & l> 0.15)
stem <-  makeStem(ht.max,dbh/100,0.05,15)
cstem <- stem
for(i in 1:nrow(shapes)){#i=1
  branch <- makeStem(shapes$l[i], shapes$d[i]*.5,0.01,10)
  bpos <- max(branch$y)
  branch2 <- branch |> mutate(x=x*0.3*bf,y=y*0.3*bf)
  branch3 <- branch |> mutate(x=x*0.2*bf,y=y*0.2*bf)
  branchx1 <- branch |> mutate(x=x*0.2*bf,y=y*0.2*bf)
  branch <- skewStem(branch, amp=ifelse(shapes$a[i] >= 0,-0.07*(shapes$s[i]*-1+1),0.07*(shapes$s[i]*-1+1)),
                     phase=0, waves=1)
  branchA <- attachBranch(branch, branch2, ifelse(shapes$a[i] >= 0,30,-30), bpos*0.3)#lower
  branchA <- attachBranch(branchA, branch3, ifelse(shapes$a[i] >= 0,-30,30), bpos*0.7)#upper
  #basal branch twigs to attach crown
  branchB <- attachBranch(branchA, branchx1, ifelse(shapes$a[i] >= 0,-90,90), bpos*0.05)
  branchB <- attachBranch(branchB, branchx1, ifelse(shapes$a[i] >= 0,50,-50), bpos*0.05)
  
  stem <- attachBranch(stem, branchA, shapes$a[i], shapes$by[i])#branches to show bare
  cstem <- attachBranch(cstem, branchB, shapes$a[i], shapes$by[i])#branches to attach crown
}
crown <- cstem |> subset(grepl('tip',type))
crown2 <- vegnasis::cavhull(x=crown$x,y=crown$y, concave = F)

 ggplot()+
  geom_polygon(data=stem, aes(x=x, y=y), color='brown',fill='#99500090')+
  # geom_point(data=stem, aes(x=x, y=y), color='red')+
   geom_polygon(data=crown2, aes(x=x, y=y), color='green',fill='#00990090')+
   geom_polygon(data=crown, aes(x=x, y=y), color='green',fill='#00990090')+
  # geom_point(data=tree, aes(x=x, y=y), color='green')+
  coord_fixed()
