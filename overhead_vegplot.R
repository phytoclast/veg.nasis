library(vegnasis)
library(ggplot2)

veg.raw <-  vegnasis::nasis.veg
veg0 <- clean.veg(veg.raw)

veg <- subset(veg0,  grepl('2022MI165021.P',plot))

veg <- veg  |> fill.type.df() |> fill.hts.df()
#check to see if habit code exists, then only fill in missing values---
if(!'habit'%in% colnames(veg)){veg <- veg |> mutate(habit = get.habit.code(taxon))
}else{#override default values
  veg <- veg |> mutate(prehabit =get.habit.code(taxon),
                       habit = ifelse(is.na(veg$habit),prehabit,habit),
                       prehabit=NULL)}
veg <- veg |> mutate(dbh.r =  fill.diameters(ht.max,dbh.max,dbh.min))
veg <- veg |> mutate(cw =  case_when(grepl('^T', habit) | ht.max > 5 ~ pmax(est_crown_width(dbh.r),1),
                                     grepl('^S', habit) ~ pmax(pmin(3,ht.max),1),
                                     TRUE ~ 1))
veg <- veg |> mutate(density =  density_from_cw(cover, cw))
veg <- veg |> mutate(BA.r =  BA_per_ha(density, dbh.r))

veg <- veg |> group_by(plot) |> mutate(BA.sum = sum(BA, na.rm = T), BA.rsum = sum(BA.r, na.rm = T), BA.sum = ifelse(is.na(BA.sum) | BA.sum <=0, BA.rsum,BA.sum), BA.ratio = BA.sum/BA.rsum,  BA.rsum = NULL)#

veg <- veg |> mutate(BA.r = ifelse(ht.max <= 5, BA.r, round(BA.r*BA.ratio,1)),
                     density = ifelse(ht.max <= 5, density, round(density*BA.ratio,0)),
                     cw = ifelse(ht.max <= 5, cw, round(cw*BA.ratio^-0.5,1)))

veg <- veg |> mutate(crshape0 = case_when(grepl('^T', habit) & grepl('NE', habit) ~ 'conifer1',
                                          grepl('^T', habit) & grepl('N', habit) ~ 'conifer3',
                                          grepl('^T', habit) & grepl('P', habit) ~ 'palm',
                                          grepl('^T', habit) & grepl('F', habit) ~ 'palm',
                                          grepl('^T', habit) & grepl('BE', habit) ~ 'blob2',
                                          grepl('^T', habit)  ~ 'blob1',
                                          grepl('^S', habit) & grepl('P', habit) ~ 'palm',
                                          grepl('^S', habit) & grepl('F', habit) ~ 'palm',
                                          grepl('^S', habit) & grepl('BE', habit) ~ 'blob2',
                                          grepl('^S', habit) ~ 'cloud1',
                                          grepl('FE', habit) ~ 'ferny',
                                          grepl('^H', habit) & grepl('F', habit) ~ 'forby',
                                          grepl('^H', habit) & type %in% 'grass/grasslike' ~ 'grassy'),
                     crfill0 = case_when(grepl('^T', habit) & grepl('NE', habit) ~ '#1A801A',
                                         grepl('^T', habit) & grepl('N', habit) ~ 'green',
                                         grepl('^T', habit) & grepl('P', habit) ~ 'green',
                                         grepl('^T', habit) & grepl('F', habit) ~ 'green',
                                         grepl('^T', habit) & grepl('BE', habit) ~ '#1A801A',
                                         grepl('^T', habit)  ~ '#4DE600',
                                         grepl('^S', habit) & grepl('P', habit) ~ 'green',
                                         grepl('^S', habit) & grepl('F', habit) ~ 'green',
                                         grepl('^S', habit) & grepl('BE', habit) ~ '#1A801A',
                                         grepl('^S', habit) ~ 'green',
                                         grepl('FE', habit) ~ 'green',
                                         grepl('^H', habit) & grepl('F', habit) ~ 'magenta',
                                         grepl('^H', habit) & type %in% 'grass/grasslike' ~ 'yellowgreen'),
                     crcolor0 = case_when(grepl('^T', habit) & grepl('NE', habit) ~ 'darkgreen',
                                          grepl('^T', habit) & grepl('N', habit) ~ 'darkgreen',
                                          grepl('^T', habit) & grepl('P', habit) ~ 'darkgreen',
                                          grepl('^T', habit) & grepl('F', habit) ~ 'darkgreen',
                                          grepl('^T', habit) & grepl('BE', habit) ~ 'darkgreen',
                                          grepl('^T', habit)  ~ 'darkgreen',
                                          grepl('^S', habit) & grepl('P', habit) ~ 'darkgreen',
                                          grepl('^S', habit) & grepl('F', habit) ~ 'darkgreen',
                                          grepl('^S', habit) & grepl('BE', habit) ~ 'darkgreen',
                                          grepl('^S', habit) ~ 'darkgreen',
                                          grepl('FE', habit) ~ 'darkgreen',
                                          grepl('^H', habit) & grepl('F', habit) ~ 'darkgreen',
                                          grepl('^H', habit) & type %in% 'grass/grasslike' ~ '#4D8000'),
                     
                     stshape0 = case_when(grepl('^T', habit) & grepl('N', habit) ~ 'trunk',
                                          grepl('^T', habit)  ~ 'trunk',
                                          grepl('^S', habit)  ~ 'sticks',
                                          grepl('FE', habit) ~ NA,
                                          grepl('F', habit) ~ NA,
                                          type %in% 'grass/grasslike' ~ NA),
                     stfill0 = case_when(grepl('^T', habit) & grepl('N', habit) ~ 'orange',
                                         grepl('^T', habit)  ~ 'orange',
                                         grepl('^S', habit)  ~ 'orange',
                                         grepl('FE', habit) ~ NA,
                                         grepl('F', habit) ~ NA,
                                         type %in% 'grass/grasslike' ~ NA),
                     stcolor0 = case_when(grepl('^T', habit) & grepl('N', habit) ~ 'brown',
                                          grepl('^T', habit)  ~ 'brown',
                                          grepl('^S', habit)  ~ 'brown',
                                          grepl('FE', habit) ~ NA,
                                          grepl('F', habit) ~ NA,
                                          type %in% 'grass/grasslike' ~ NA),
                     fun=case_when(grepl('^T', habit)  ~ 'T',
                                   grepl('^S', habit)  ~ 'S',
                                   grepl('^H', habit)  ~ 'H'),
                     stems = round(plength*pwidth/10000*density,0) #count number of stem objects required plot size.
)
#Check if user defined values exist for shape, fill, and color.
if('crshape' %in% colnames(veg)){
  veg <- veg |> mutate(crshape = ifelse(is.na(crshape), crshape0, crshape), crshape0=NULL)}else{
    veg <- veg |> mutate(crshape = crshape0, crshape0=NULL)}
if('crfill' %in% colnames(veg)){
  veg <- veg |> mutate(crfill = ifelse(is.na(crfill), crfill0, crfill), crfill0=NULL)}else{
    veg <- veg |> mutate(crfill = crfill0, crfill0=NULL)}
if('crcolor' %in% colnames(veg)){
  veg <- veg |> mutate(crcolor = ifelse(is.na(crcolor), crcolor0, crcolor), crcolor0=NULL)}else{
    veg <- veg |> mutate(crcolor = crcolor0, crcolor0=NULL)}

if('stshape' %in% colnames(veg)){
  veg <- veg |> mutate(stshape = ifelse(is.na(stshape), stshape0, stshape), stshape0=NULL)}else{
    veg <- veg |> mutate(stshape = stshape0, stshape0=NULL)}
if('stfill' %in% colnames(veg)){
  veg <- veg |> mutate(stfill = ifelse(is.na(stfill), stfill0, stfill), stfill0=NULL)}else{
    veg <- veg |> mutate(stfill = stfill0, stfill0=NULL)}
if('stcolor' %in% colnames(veg)){
  veg <- veg |> mutate(stcolor = ifelse(is.na(stcolor), stcolor0, stcolor), stcolor0=NULL)}else{
    veg <- veg |> mutate(stcolor = stcolor0, stcolor0=NULL)}


strats <- veg |> subset(fun %in% c('T','S','H') & stems > 0) |> arrange(plot,-ht.max, -cover)#reduce to strata which have models and are not empty of stems
strats$seq <- c(1:nrow(strats))
strats <- strats |> group_by(plot) |> mutate(seqmin = min(seq), seq = seq-seqmin+1, seqmin = NULL)


#Establish a stand 50 by 20 m.
stand <- make_hex_stand(plength/100,1) |> subset(yp >= mean(yp)-pwidth/2 & yp < mean(yp)+pwidth/2) |> mutate(wtn = wt, stratid = NA)

#assign stump locations per stratum
for (i in 1:nrow(strats)){#i=1
  thistrat = strats$seq[i]
  nstems = strats$stems[i]
  newstumps <- sample(stand$stumpid, size = nstems, prob = stand$wtn, replace = T)
  stand <- stand |> mutate(wtn = ifelse(stand$stumpid %in% newstumps, 0.0001, wtn),
                           stratid = ifelse(stand$stumpid %in% newstumps, thistrat,stratid))
}





angles <- c(0:59)/60*2*3.141592
#coniferoverhead
wave = ((sin(angles*9)+1)/2)^1
x = cos(angles)*(wave+1)/2; y=sin(angles)*(wave+1)/2
shape='conifer'
conifer <- data.frame(shape=shape, x=x,y=y)
#hardwoodoverhead
wave = ((sin(angles*9)+1)/2)^0.5
x = cos(angles)*(wave+4)/5; y=sin(angles)*(wave+4)/5
shape='hardwood'
hardwood <- data.frame(shape=shape, x=x,y=y)
#circle
x = cos(angles); y=sin(angles)
shape='circle'
circle <- data.frame(shape=shape, x=x,y=y)
overheadshapes <- rbind(circle, hardwood, conifer)



ggplot(overheadshapes, aes(x=x, y=y, color=shape))+
  geom_polygon(fill='#FFFFFF00')


make_tree.overhead <- function(ht.max, ht.min, crwd, dbh, crshape, stshape){
  crown <- subset(overheadshapes, shape %in% 'hardwood') |> mutate(x=x*crwd, y=y*crwd, obj='crown', ht.max=ht.max)
  base <- subset(overheadshapes, shape %in% 'circle') |> mutate(x=x*dbh/100*1.1, y=y*dbh/100*1.1, obj='stem', ht.max=ht.max)
  tree = rbind(crown, base)
  tree$ptord <- rownames(tree) |> as.numeric()
  return(tree)}

make_plant_overhead<- function(fun, ht.max, ht.min,crwd,dbh, crshape, stshape){

    plant <- make_tree.overhead(ht.max, ht.min, crwd, dbh, crshape, stshape)
  
  plant <- plant |> mutate(fill=NULL,color=NULL)
  return(plant)
}









#Create shapes of the right size, then distribute into the stump positions.
for (i in 1:nrow(strats)){#i=1
  thistrat <- strats[i,]
  plant0 <- make_plant_overhead(thistrat$fun, thistrat$ht.max, thistrat$ht.min,thistrat$cw,thistrat$dbh.r, thistrat$crshape, thistrat$stshape)
  stumps0 <- stand |> subset(stratid %in% i)
  plant0 <- merge(stumps0, plant0) |> mutate(objid = paste0(stratid,obj,stumpid))
  colors0 <- subset(thistrat, select=c(crfill, stfill, crcolor, stcolor))
  plant0 <- merge(plant0,colors0) |> mutate(fill = ifelse(obj %in% c('crown','herb'), crfill, stfill),
                                            color = ifelse(obj %in% c('crown','herb'), crcolor, stcolor),
                                            crfill = NULL, stfill = NULL, crcolor = NULL, stcolor = NULL)
  if(i==1){plants <- plant0}else{plants <- rbind(plants,plant0)}
}

#randomize sizes and positions
plants <- plants |> group_by(stumpid) |>
  mutate(crwd = max(x)-min(x),
         xpp = xp + runif(1, min = -0.8, max = 0.8),#shift position on grid
         ypp = yp + runif(1, min = -0.8, max = 0.8),#shift position on grid
         xr = rnorm(1,crwd, crwd/10)/crwd,#deviation in width
         xn = x*xr+xpp,#resized width and put on new position
         yn = y*xr+ypp)





















plength = 50; pwidth=50
ytrans = 'identity'; yratio=1; units = 'm'; skycolor = "#D9F2FF80"; fadecolor = "#D9F2FF"; gridalpha=0.3; groundcolor="#808066"; xlim=c(0,50); ylim=c(0,50); xticks=5; yticks=5; xslope=0; yslope=0; xperiod=10; xamplitude=0


xnmax <- max(plants$xn, na.rm =TRUE)
xnmin <- min(plants$xn, na.rm =TRUE)
ynmax <- max(plants$yn, na.rm =TRUE)
ynmin <- min(plants$yn, na.rm =TRUE)

iters = round(max(plants$ht.max)+1,0) #total veg height
groundline = data.frame(xn=c(xnmin,xnmax,xnmax,xnmin),
                        yn=c(ynmin,ynmin,ynmax,ynmax))
                        
ground = data.frame(groundline)
pcolor <- c(plants$color, ground$color) |> unique() |> sort()
pfill <- c(plants$fill, ground$fill) |> unique()|> sort()

ucf = case_when(units %in% c('feet', 'ft') ~ 0.3048,
                units %in% c('inches', 'in') ~ 0.3048/12,
                units %in% c('cm') ~ 0.01,
                TRUE ~ 1)
units = ifelse(ucf == 1, 'm',units)

yunits = paste0('ground distance (', units,')')
xunits = paste0('ground distance (', units,')')
ybreaks = seq(floor(ylim[1]/ucf/yticks)*yticks-yticks,floor(ylim[2]/ucf/yticks)*xticks+yticks,yticks)*ucf
xbreaks = seq(floor(xlim[1]/ucf/xticks)*xticks-xticks,floor(xlim[2]/ucf/xticks)*xticks+xticks,xticks)*ucf
yminor = seq(floor(ylim[1]/ucf-yticks),floor(ylim[2]/ucf+yticks),yticks/5)*ucf
xminor = seq(floor(xlim[1]/ucf-xticks),floor(xlim[2]/ucf+xticks),xticks/5)*ucf
ylabels = ybreaks/ucf
xlabels =  xbreaks/ucf





gp <- ggplot()+
  geom_polygon(data=ground, aes(x=xn,y=yn), fill=groundcolor, color=groundcolor, alpha=1, linewidth=0.1)


  for(i in 1:iters){#i=25 i=1:iters
  plants0 <- plants |> subset(ht.max <= i &  ht.max > i-1)
  if(nrow(plants0)>0){stems0 <- plants0 |> subset(obj %in% 'stem')
  crowns0 <- plants0 |> subset(obj %in% c('crown','herb'))
  gp = gp+
    
    geom_polygon(data=crowns0, aes(x=xn,y=yn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
    geom_polygon(data=stems0, aes(x=xn,y=yn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)}
}
gp = gp +
  scale_fill_manual(values=pfill)+
  scale_color_manual(values=pcolor)+
  theme(legend.position = "none",
        
        panel.background = element_rect(fill = skycolor,
                                        colour = "black",
                                        linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                        colour = rgb(0.1, 0.1, 0.1, gridalpha)),
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid',
                                        colour = rgb(0.1, 0.1, 0.1, gridalpha/3))
  )+
  coord_fixed(ratio = yratio, ylim=ylim,xlim=xlim, expand = FALSE)+
  scale_y_continuous(name = yunits, trans = ytrans, labels = ylabels, breaks = ybreaks, minor_breaks = yminor, limits = c(ynmin-5,ynmax+5))+#
  scale_x_continuous(name = xunits ,breaks = xbreaks, labels = xlabels, minor_breaks = xminor, limits = c(xnmin-5,xnmax+5))#
