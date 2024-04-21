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
                         d1 = floor(d)) |> group_by(d1, y>=0) |> 
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