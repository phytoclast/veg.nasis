library(soilDB)
library(vegnasis)
library(cluster)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

siteass <- get_site_association_from_NASIS(SS=F)
sites <- get_site_data_from_NASIS_db(SS=F)
veg.raw <- soilDB::get_vegplot_species_from_NASIS_db(SS=F)
vegplot <- soilDB::get_vegplot_from_NASIS_db(SS=F)
vegground <- get_vegplot_groundsurface_from_NASIS_db(SS=F)
ba <- get_vegplot_speciesbasalarea_from_NASIS(SS=F)
si <- get_vegplot_tree_si_details_from_NASIS_db(SS=F)
trans <- get_vegplot_transect_from_NASIS_db(SS=F)
transp <- get_vegplot_transpecies_from_NASIS_db(SS=F)

ecosites <- sites[,c('usiteid','ecositeid','ecositenm','commphasename')] |> left_join(vegplot[,c('usiteid','vegplotid')], relationship = 'many-to-many') |> subset(!is.na(vegplotid) & !is.na(ecositeid)) |> unique()

veg0 <- veg.raw |> clean.veg() |> subset(!is.na(plot))

#select ecosites ----
ecosites <- sites[,c('usiteid','ecositeid','ecositenm','commphasename')] |> left_join(vegplot[,c('usiteid','vegplotid')], relationship = 'many-to-many') |> subset(!is.na(vegplotid) & !is.na(ecositeid)) |> unique()



#select ecosites ----
thoseeco <- c('F097XA006MI','F097XA007MI','F096XB019MI','F096XB020MI','F096XB021MI','F096XB023MI','F096XA007MI','F096XA008MI')
# thoseeco <- c('F096XB020MI')

thiseco <- ecosites |> subset(ecositeid %in% thoseeco)

veg <- veg0 |> subset(plot %in% thiseco$vegplotid)

#Create alternative higher taxon datasets and combine them.
veg.genera <- veg |> mutate(taxon = link.taxonomy(taxon, taxrank=1))
veg.families <- veg |> mutate(taxon = link.taxonomy(taxon, taxrank=2))
veg.combined <- rbind(veg,veg.genera,veg.families)

#Create plot matrix, based log transformed relative cover values.
m <- veg.combined |> make.plot.matrix(tr = 'log', rc = TRUE)




#distance matrix based on Bray-Curtis simularity.
d = vegan::vegdist(m, method='bray')
#Cluster analysis using Ward's method using distance matrix.
t <- cluster::agnes(d, method = 'gaverage', par.method = -0.10)|> as.hclust()
#Define number of groups to color the dendrogram by.
k = 10
groups <- cutree(t, k = k)
#This function rearranges the branchs and groups so that the tree is always oriented with most nested branches to the bottom of the plot (when tree oriented vertically with branches to the right).
groups <- dendrogrouporder(t, groups)
a = 'Vegetation'

plot.dendro(a,d,t,groups)
######find optimal clusters ----
meths <- c('nmds_pam', 'nmds_kmeans', 'bray_pam', 'bray_kmeans', 'bray_ward', 'bray_flex','nmds_flex')
ks <- 2:15
betas <- seq(from = -1, to = 1, length.out = 41) |> round(2)
df <- NULL
for(i.m in 1:length(meths)){#i.m=1
  thismeth <-  meths[i.m]
  
  if(thismeth %in% 'bray_ward'){
    t <- cluster::agnes(d, method = 'ward')|> as.hclust()
    for(i.k in 1:length(ks)){
      k=ks[i.k]
      c1 <- cutree(t, k = k)
      sil <- cluster::silhouette(c1, dist = d)
      silmean <- mean(sil[,3])
      silmed <- median(sil[,3])
      silmin <- quantile(sil[,3],0.25)
      df0 <- data.frame(meth = thismeth, k = k, bet = NA, silmean = silmean, silmed = silmed, silmin = silmin)
      if(is.null(df)){df <- df0}else{df <- rbind(df,df0)}
    }
  }else if(thismeth %in% 'bray_flex'){
    for(i.b in 1:length(betas)){#i.b=1
      bet <- betas[i.b]
      t <- cluster::agnes(d, method = 'gaverage', par.method = bet)|> as.hclust()
      for(i.k in 1:length(ks)){
        k=ks[i.k]
        c1 <- cutree(t, k = k)
        sil <- cluster::silhouette(c1, dist = d)
        silmean <- mean(sil[,3])
        silmed <- median(sil[,3])
        silmin <- quantile(sil[,3],0.25)
        df0 <- data.frame(meth = thismeth, k = k, bet = bet, silmean = silmean, silmed = silmed, silmin = silmin)
        if(is.null(df)){df <- df0}else{df <- rbind(df,df0)}
      }}
  }else if(thismeth %in% 'nmds_flex'){
    for(i.b in 1:length(betas)){#i.b=1
      bet <- betas[i.b]
      t <- cluster::agnes(mtx, method = 'gaverage', par.method = bet)|> as.hclust()
      for(i.k in 1:length(ks)){
        k=ks[i.k]
        c1 <- cutree(t, k = k)
        sil <- cluster::silhouette(c1, dist = d)
        silmean <- mean(sil[,3])
        silmed <- median(sil[,3])
        silmin <- quantile(sil[,3],0.25)
        df0 <- data.frame(meth = thismeth, k = k, bet = bet, silmean = silmean, silmed = silmed, silmin = silmin)
        if(is.null(df)){df <- df0}else{df <- rbind(df,df0)}
      }}
  }else if(thismeth %in% 'nmds_pam'){
    for(i.k in 1:length(ks)){
      k=ks[i.k]
      c1 <- mtx |> pam(k, nstart = 20) 
      c1 <- c1$cluster
      sil <- cluster::silhouette(c1, dist = d)
      silmean <- mean(sil[,3])
      silmed <- median(sil[,3])
      silmin <- quantile(sil[,3],0.25)
      df0 <- data.frame(meth = thismeth, k = k, bet = NA, silmean = silmean, silmed = silmed, silmin = silmin)
      if(is.null(df)){df <- df0}else{df <- rbind(df,df0)}
    }
  }else if(thismeth %in% 'nmds_kmeans'){
    for(i.k in 1:length(ks)){
      k=ks[i.k]
      c1 <- mtx |> pam(k, nstart = 20) 
      c1 <- c1$cluster
      sil <- cluster::silhouette(c1, dist = d)
      silmean <- mean(sil[,3])
      silmed <- median(sil[,3])
      silmin <- quantile(sil[,3],0.25)
      df0 <- data.frame(meth = thismeth, k = k, bet = NA, silmean = silmean, silmed = silmed, silmin = silmin)
      if(is.null(df)){df <- df0}else{df <- rbind(df,df0)}
    }
  }else if(thismeth %in% 'bray_kmeans'){
    for(i.k in 1:length(ks)){
      k=ks[i.k]
      c1 <- d |> kmeans(k, nstart = 20) 
      c1 <- c1$cluster
      sil <- cluster::silhouette(c1, dist = d)
      silmean <- mean(sil[,3])
      silmed <- median(sil[,3])
      silmin <- quantile(sil[,3],0.25)
      df0 <- data.frame(meth = thismeth, k = k, bet = NA, silmean = silmean, silmed = silmed, silmin = silmin)
      if(is.null(df)){df <- df0}else{df <- rbind(df,df0)}
    }
  }else if(thismeth %in% 'bray_pam'){
    for(i.k in 1:length(ks)){
      k=ks[i.k]
      c1 <- d |> pam(k, nstart = 20)
      c1 <- c1$cluster
      sil <- cluster::silhouette(c1, dist = d)
      silmean <- mean(sil[,3])
      silmed <- median(sil[,3])
      silmin <- quantile(sil[,3],0.25)
      df0 <- data.frame(meth = thismeth, k = k, bet = NA, silmean = silmean, silmed = silmed, silmin = silmin)
      if(is.null(df)){df <- df0}else{df <- rbind(df,df0)}
    }}}

# ggplot(df, aes(x=as.factor(k),y=silmean))+
#   geom_boxplot()
# ggplot(df, aes(x=meth,y=silmin))+
#   geom_boxplot()
# ggplot(subset(df, meth %in% 'bray_flex'), aes(x=as.factor(bet),y=silmean))+
#   geom_boxplot()
# ggplot(subset(df, meth %in% 'bray_flex'), aes(x=as.factor(bet),y=silmed))+
#   geom_boxplot()
# ggplot(subset(df, meth %in% 'nmds_flex'), aes(x=as.factor(bet),y=silmin))+
#   geom_boxplot()

# ggplot(subset(df, meth %in% c('bray_flex','nmds_flex') & bet == -0.1 | meth %in% c('bray_kmeans','bray_ward','nmds_kmeans','bray_pam','nmds_pam')), aes(x=k,y=silmin, color = meth))+
#   geom_line(linewidth = 1)+
#   geom_point(size= 2)+
#   scale_x_continuous(breaks = 2:15)

bestbet <- subset(df, meth %in% 'bray_flex') |> group_by(bet) |> summarise(sil = mean(silmin))
ggplot(subset(df, meth %in% c('bray_flex') & bet == -0.10 | meth %in% c('bray_kmeans','bray_ward','bray_pam')), aes(x=k,y=silmed, color = meth))+
  geom_line(linewidth = 1)+
  geom_point(size= 2)+
  scale_x_continuous(breaks = 2:15)

# ggplot(subset(df, meth %in% c('nmds_flex') & bet == -.1 | meth %in% c('nmds_kmeans','bray_ward','nmds_pam')), aes(x=k,y=silmin, color = meth))+
#   geom_line(linewidth = 1)+
#   geom_point(size= 2)+
#   scale_x_continuous(breaks = 2:15)


#Nonmetric Multidimensional Scaling (NMDS) ----
library(vegan)
ndim <- 3 #number of dimensions
set.seed(6190)
nmds <- vegan::metaMDS(m, k=ndim, try=50, trymax = 50)

#join NMDS data components
spts.es <- data.frame(plot=rownames(m)) |> cbind(m)  
pt.df <- vegan::scores(nmds, display='sites') |> as_tibble(rownames='sites')  
sp.df <- vegan::scores(nmds, display='species') |> as_tibble(rownames='species')

#narrow list of species and assign PLANTS Symbols for compact display on plot.
spsums <- veg |> group_by(plot, taxon) |> summarise(cover = vegnasis::cover.agg(cover)) |> subset(cover > 1) |> group_by(taxon) |> summarise(nplot = length(cover), maxcover = max(cover)) |> subset(nplot >= 2 | maxcover >=25)
sp.df <- sp.df |> mutate(sp1 = stringr::str_replace_all(species,'\\.', ' '), sp2 = vegnasis::fill.usda.symbols(sp1)) |> subset(sp1 %in% spsums$taxon)

#cut into clusters
mtx <- scores(nmds, display='sites')
k <- 10
set.seed(6190)
pam1 <- d |> pam(k, nstart = 20)
pam1<-pam1$clustering
tf <- cluster::agnes(d, method = 'gaverage', par.method = -0.1)|> as.hclust()
flex1 <- cutree(tf, k = k)
flex1 <- dendrogrouporder(tf, flex1)
tw <- cluster::agnes(d, method = 'ward')|> as.hclust()
ward1 <- cutree(tw, k = k)
ward1 <- dendrogrouporder(tw, ward1)
# plot.dendro('ward',d,tw,ward1)
pt.df <- pt.df |> mutate(pamclust = as.factor(pam1), wardclust = as.factor(ward1), flexclust = as.factor(flex1)) 

#convex hulls
fromclusts <- unique(pt.df$pamclust)
for(i in 1:k){#i=1
  thisclust <-fromclusts[i]
  thishull <- pt.df |> subset(pamclust %in% thisclust)
  chull0 <- chull(thishull$NMDS1, thishull$NMDS2)
  thishull  <- thishull[chull0,]
  if(i==1){pamhull=thishull}else{pamhull=rbind(pamhull,thishull)}
};rm(thisclust,thishull,chull0)

fromclusts <- unique(pt.df$wardclust)
for(i in 1:k){#i=1
  thisclust <-fromclusts[i]
  thishull <- pt.df |> subset(wardclust %in% thisclust)
  chull0 <- chull(thishull$NMDS1, thishull$NMDS2)
  thishull  <- thishull[chull0,]
  if(i==1){wardhull=thishull}else{wardhull=rbind(wardhull,thishull)}
};rm(thisclust,thishull,chull0)

fromclusts <- unique(pt.df$flexclust)
for(i in 1:k){#i=1
  thisclust <-fromclusts[i]
  thishull <- pt.df |> subset(flexclust %in% thisclust)
  chull0 <- chull(thishull$NMDS1, thishull$NMDS2)
  thishull  <- thishull[chull0,]
  if(i==1){flexhull=thishull}else{flexhull=rbind(flexhull,thishull)}
};rm(thisclust,thishull,chull0)
######
fromclusts <- unique(pt.df$pamclust)
for(i in 1:k){#i=1
  thisclust <-fromclusts[i]
  thishull <- pt.df |> subset(pamclust %in% thisclust)
  chull0 <- chull(thishull$NMDS1, thishull$NMDS3)
  thishull  <- thishull[chull0,]
  if(i==1){pamhull2=thishull}else{pamhull2=rbind(pamhull2,thishull)}
};rm(thisclust,thishull,chull0)

fromclusts <- unique(pt.df$wardclust)
for(i in 1:k){#i=1
  thisclust <-fromclusts[i]
  thishull <- pt.df |> subset(wardclust %in% thisclust)
  chull0 <- chull(thishull$NMDS1, thishull$NMDS3)
  thishull  <- thishull[chull0,]
  if(i==1){wardhull2=thishull}else{wardhull2=rbind(wardhull2,thishull)}
};rm(thisclust,thishull,chull0)

fromclusts <- unique(pt.df$flexclust)
for(i in 1:k){#i=1
  thisclust <-fromclusts[i]
  thishull <- pt.df |> subset(flexclust %in% thisclust)
  chull0 <- chull(thishull$NMDS1, thishull$NMDS3)
  thishull  <- thishull[chull0,]
  if(i==1){flexhull2=thishull}else{flexhull2=rbind(flexhull2,thishull)}
};rm(thisclust,thishull,chull0)


library(ggplot2)

gp <- ggplot() +
  # geom_polygon(data=pamhull, aes(x=NMDS1,y=NMDS2, color=pamclust), alpha=0, linewidth=1)+
  # geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2, color=pamclust), alpha=0.5, size=2)+
  geom_polygon(data=flexhull, aes(x=NMDS1,y=NMDS2, color=flexclust), alpha=0, linewidth=1)+
  geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2, color=flexclust), alpha=0.5, size=2)+
  # geom_polygon(data=wardhull, aes(x=NMDS1,y=NMDS2, color=wardclust), alpha=0, linewidth=1)+
  # geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2, color=wardclust), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS2), color='blue')+
  geom_text(data=sp.df, aes(label=sp2, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue', size=2)+
  coord_fixed()
gp

gp <- ggplot() +
  # geom_polygon(data=pamhull2, aes(x=NMDS1,y=NMDS3, color=pamclust), alpha=0, linewidth=1)+
  # geom_point(data=pt.df, aes(x=NMDS1,y=NMDS3, color=pamclust), alpha=0.5, size=2)+
  geom_polygon(data=flexhull2, aes(x=NMDS1,y=NMDS3, color=flexclust), alpha=0, linewidth=1)+
  geom_point(data=pt.df, aes(x=NMDS1,y=NMDS3, color=flexclust), alpha=0.5, size=2)+
  # geom_polygon(data=wardhull2, aes(x=NMDS1,y=NMDS3, color=wardclust), alpha=0, linewidth=1)+
  # geom_point(data=pt.df, aes(x=NMDS1,y=NMDS3, color=wardclust), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS3), color='blue')+
  geom_text(data=sp.df, aes(label=sp2, x=NMDS1,y=NMDS3), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue', size=2)+
  coord_fixed()
gp




  
  
  #Describe clusters ----

  veg <- veg |> fill.type.df() |> fill.hts.df()
  association <- veg |> get.assoc()
  struct <- veg |> get.structure(simple = FALSE)
  # k=k
  # set.seed(6190)
  # pam1 <- d |> pam(k, nstart = 20)
  # pam1<-pam1$clustering
  # tf <- cluster::agnes(d, method = 'gaverage', par.method = -0.1)|> as.hclust()
  # flex1 <- cutree(tf, k = k)
  # flex1 <- dendrogrouporder(tf, flex1)
  # tw <- cluster::agnes(d, method = 'ward')|> as.hclust()
  # ward1 <- cutree(tw, k = k)
  # ward1 <- dendrogrouporder(tw, ward1)
  mean(cluster::silhouette(pam1, dist = d)[,3])
  mean(cluster::silhouette(flex1, dist = d)[,3])
  mean(cluster::silhouette(ward1, dist = d)[,3])
  cldf <- data.frame(plot=names(pam1), pam = (pam1))
  flexdf <- data.frame(plot=names(flex1), flex = (flex1))
  warddf <- data.frame(plot=names(ward1), ward = (ward1))
  
  association1 <- association |> left_join(struct) |> left_join(ecosites[,c("ecositeid","ecositenm", "vegplotid" , "commphasename")], join_by(plot==vegplotid)) |> left_join(cldf) |> left_join(warddf) |> left_join(flexdf) 
  
  
write.csv(association1, 'all-F096XB021MI.csv', na='', row.names = F)

  plot.dendro('flex-beta -0.1',d,tf,flex1)
  plot.dendro('ward',d,tw,ward1)
  

  
  library(sf)
  library(mapview)
  myplots <- veg |> left_join(ecosites[,c("usiteid","vegplotid")], join_by(plot == vegplotid))
  s <- sites |> subset(usiteid %in% myplots$usiteid) |>
    mutate(lat=latstddecimaldegrees, lon = longstddecimaldegrees) |> subset(!is.na(lon), select=c(site_id, obsdate, lat, lon, elev, slope, aspect, site_mlra, site_state, site_county, ecositeid, ecositenm, ecostatename, commphasename)) |> mutate(commphasename = tolower(commphasename)) |> 
    subset(site_id %in% 'S2025MI101001')
    #subset(commphasename %in% tolower(c('Transitional Pine Plantation')))
  s <- s |> st_as_sf(coords = c(x='lon', y='lat'), crs=st_crs('EPSG:4326'))
  
  #col.regions=c('red', 'green', 'yellow','orange','blue', 'violet', 'white','black')
  mapview(s, col.regions=c('red', 'green', 'yellow','orange','blue', 'violet', 'white','black'), zcol='site_id')
  
  # tolower(c('Oak-Pine Barrens','Oak Barrens','Native Ruderal Forest','Red Pine Plantation','Transitional Pine Plantation','White Pine Plantation','Xerophytic Forest','Disturbed Forest','Native Ruderal Forest'))
  