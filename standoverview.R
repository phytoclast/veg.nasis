library(ggplot2)
library(vegnasis)

shapes <- vegnasis::shapes

#Establish a stand 50 by 20 m.
plength=100; pwidth=100
stand <- make_hex_stand(plength/100,1) |> subset(yp >= mean(yp)-pwidth/2 & yp < mean(yp)+pwidth/2) |> mutate(wtn = wt, stratid = NA)

veg.raw <-  vegnasis::nasis.veg
veg <- clean.veg(veg.raw)
veg0 <- veg  |> fill.type.df() |> fill.hts.df()
vegstr <-  get.structure(veg0)
vegstr <- vegstr |> mutate(pcover = NA)
isforest <- subset(vegstr, tree > 0)
veg0 <- subset(veg0,  plot %in% isforest$plot)
projs <- veg0$plot |> unique()
for (j in 1:length(projs)) {#j=1
  thisproj <- projs[j]
veg <- subset(veg0,  plot %in% thisproj)







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


#assign stump locations per stratum
for (i in 1:nrow(strats)){#i=1
  thistrat = strats$seq[i]
  nstems = strats$stems[i]
  newstumps <- sample(stand$stumpid, size = nstems, prob = stand$wtn, replace = T)
  stand <- stand |> mutate(wtn = ifelse(stand$stumpid %in% newstumps, 0.0001, wtn),
                           stratid = ifelse(stand$stumpid %in% newstumps, thistrat,stratid))
}



standoverview <- stand |> left_join(strats[,c('crshape','dbh.r','cw','ht.max','seq')], by=join_by(stratid==seq))
blob <- subset(shapes, shape %in% 'blob1') |> mutate(y = z-0.5)


strats.overstory <- standoverview |> subset(ht.max > 5)

strats.overstory <- strats.overstory |> merge(blob)

strats.overstory <- strats.overstory |> group_by(stumpid) |> mutate(nxp = xp+rnorm(1,0, 1/2), nyp = yp+rnorm(1,0, 1/2), ncw = rnorm(1,cw, cw/10), xn = nxp+x*ncw, yn=nyp+y*ncw) |> ungroup()

gp <- ggplot() +
  geom_polygon(data=strats.overstory, aes(x=xn,y=yn,group=stumpid), fill='darkgreen', color='darkgreen', alpha=1, linewidth=0.01)+
  theme(legend.position = "none",
        panel.background = element_rect(fill = 'white',linewidth = 0),
        panel.grid.major = element_line(linewidth = 0),
        panel.grid.minor = element_line(linewidth = 0)
  )+
  coord_fixed(ylim=c(0,100),xlim=c(0,100), expand = FALSE)+
  scale_y_continuous(name = NULL, labels = NULL, breaks = NULL, minor_breaks = NULL)+
  scale_x_continuous(name = NULL, labels = NULL, breaks = NULL, minor_breaks = NULL)

library(terra)
tiff('standplot.tif', width = 1000, height = 1000)
gp
dev.off()
gprast <- rast('standplot.tif')
plot(gprast)
gprast1 <- gprast$standplot_1
gprast1 <- ifel(gprast1 <= 0,1,0)
#plot(gprast1)

gdf <- as.data.frame(gprast1)

pcover0 <- mean(gdf$standplot_1)
vegstr <- vegstr |> mutate(pcover = ifelse(plot %in% thisproj, pcover0, pcover))
}








