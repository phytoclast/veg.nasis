setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
veg.spp <- read.delim('data/Observed_Species.txt')
veg.site <- read.delim('data/Sites.txt')
veg <- clean.veg.log(veg.site,veg.spp)

veg <- veg |> fill.type.df() |> fill.hts.df()
groups <- c('2023OH051003', '2023OH051001','2023OH051002', '2023523.001', '2023523.002')
# groups <- c('s20220930.001','20230731.001', '20230801.001')
veg1 <- subset(veg, label %in% groups)



x = veg1
breaks=c(0.5,2,5,12)
woodytypes = c('tree','shrub/vine', 'epiphyte')



multi.releve.summary <-  function(x, breaks=c(0.5,5,15), woodytypes = c('tree','shrub/vine', 'epiphyte')){
  
  y <- NULL
  nbks <- length(breaks)+1
  brks <- c(0,breaks,1000)
  
  
  for(i in 1:(nbks)){#i = 5
    y0 <- x %>% subset((ht.max <= brks[i+1] & ht.max > brks[i] & type %in% woodytypes)|(i==1 & !type %in% woodytypes))
    
    if(nrow(y0)>0){
      y0 <- y0 %>% mutate(stratum=i, stratum.label = paste0(brks[i], "-", ifelse(i==nbks, "+",brks[i+1])), bottom= brks[i], top = ifelse(i==nbks,brks[i]+1,brks[i+1]))
      y1 <- y0 %>% group_by(plot, taxon, type, stratum, stratum.label, bottom, top) %>% summarise(cover = cover.agg(cover), ht.min = mean(ht.min), ht.max = mean(ht.max))
      
      if(is.null(y)){y <- y1}else{y <- rbind(y, y1)}}
  }
  
  freq0 <- y |> group_by(plot, taxon, type) |> summarise(cover = cover.agg(cover)) |> 
    mutate(freq = ifelse(cover > 0,1,0))
  nplots <- length(unique(freq0$plot))
  freq0 <- freq0 |> group_by(taxon, type) |> summarise(freq = sum(freq)/nplots)
  y <- y |> group_by(taxon, type, stratum, stratum.label, bottom, top) |> summarise(cover1 = sum(cover), ht.min = ht.round(sum(ht.min*cover)/cover1), ht.max = ht.round(sum(ht.max*cover)/cover1), cover = cover1/nplots)
  y <- y |> left_join(freq0) |> mutate(cover.p = round(cover/freq,1), cover.r = round(cover,1), freq = round(freq,3))
  y <- subset(y, select=c(taxon, type, stratum, stratum.label, freq, cover.r, cover.p, bottom, top, ht.min, ht.max))
  y <- y |> mutate(habitsort = case_when(type %in% 'tree' ~ 1,
                                         type %in% 'shrub/vine' ~ 2,
                                         type %in% 'grass/grasslike' ~ 3,
                                         type %in% 'forb' ~ 4,
                                         TRUE ~ 5))
  
  y <- y |> group_by(taxon, type) |> mutate(overstory = cover.agg(ifelse(top > 5,1,0)*cover.r), covertotal = cover.agg(cover.r)) |> arrange(habitsort, -overstory, -covertotal, taxon, -top)
  
  
  return(y)
}

x=veg1
breaks = c(0.5,2,5,12,24)
woodytypes = c('tree','shrub/vine', 'epiphyte')
lowerQ=0.25; upperQ=0.75

summary.ESIS0 <-  function(x, breaks=c(0.5,5,15), lowerQ=0.25, upperQ=0.75,woodytypes = c('tree','shrub/vine', 'epiphyte')){
  x=x |> mutate(dbh.min= ifelse(is.na(dbh.min), dbh.max,dbh.min))
  #frequency of whole plot
  f <- x |> subset(cover>0, select=c(plot,taxon)) |> unique()  |> mutate(freq=1)
  nplots = length(unique(f$plot))
  f <- f |> group_by(taxon) |> summarise(frq.plot=sum(freq)/nplots)
  
  y <- NULL
  nbks <- length(breaks)+1
  brks <- c(0,breaks,1000)
  #extract means by stratum and plot
  for(i in 1:(nbks)){#i = 1
    y0 <- x |> subset((ht.max <= brks[i+1] & ht.max > brks[i] & type %in% woodytypes)|(i==1 & !type %in% woodytypes))
    
    if(nrow(y0)>0){
      y0 <- y0 %>% mutate(stratum=i, stratum.label = paste0(brks[i], ifelse(i==nbks, "+",paste0("-", brks[i+1]))), ht.min= ht.min, ht.max = ht.max)
      y1 <- y0 %>% group_by(plot, symbol, taxon, type, stratum, stratum.label, stratum.min, stratum.max) %>%
        summarise(Cover = cover.agg(cover),
                  ht.min=weighted.mean(ht.min, cover+0.001, na.rm=TRUE),
                  ht.max=weighted.mean(ht.max, cover+0.001, na.rm=TRUE),
                  dbh.min =  weighted.mean(dbh.min, cover+0.001, na.rm=TRUE),
                  dbh.max =  weighted.mean(dbh.max, cover+0.001, na.rm=TRUE),
                  dbh.min =  ifelse(is.nan(dbh.min), NA, dbh.min),
                  dbh.max =  ifelse(is.nan(dbh.max), NA, dbh.max),
                  BA =  sum(BA, na.rm=TRUE))
      
      if(is.null(y)){y <- y1}else{y <- rbind(y, y1)}}
  }
  #weighted mean of heights of each taxon stratum among all plots
  y = y |> group_by(symbol,taxon,type,stratum, stratum.label,stratum.min, stratum.max) |>
    mutate(Bottom=round(weighted.mean(ht.min, Cover+0.001, na.rm=TRUE),1),
           Top=round(weighted.mean(ht.max, Cover+0.001, na.rm=TRUE),1),
           dbh.Low =  weighted.mean(dbh.min, Cover+0.001, na.rm=TRUE),
           dbh.High =  weighted.mean(dbh.max, Cover+0.001, na.rm=TRUE),
           dbh.Low =  ifelse(is.nan(dbh.Low), NA, round(dbh.Low,0)),
           dbh.High =  ifelse(is.nan(dbh.High), NA, round(dbh.High,0)),
           cover.ps = round(mean(Cover),1))
  
  y = y |> group_by(plot) |> mutate(totalBA = sum(BA, na.rm = TRUE), overCover = ifelse(Top > 5, Cover, NA), grossCover = sum(overCover, na.rm = TRUE), BA = round(totalBA*Cover/(grossCover+0.000001),1))
  #get frequency in stratum
  y = y |> group_by(plot,symbol,taxon,type) |> mutate(frq.strat = ifelse(sum(Cover)>0,1,0))
  #insert zeros for missing species found in other plots
  y.plot <- unique(subset(y, select=c("plot")))
  y.mid <- unique(subset(y, select=c("symbol","taxon","type","stratum","stratum.label","stratum.min", "stratum.max","Bottom","Top","dbh.Low","dbh.High","cover.ps")))
  y.fill <- merge(y.plot, y.mid) |> mutate(Cover = 0, BA = 0, frq.strat=0)
  y.fill <- rbind(y, y.fill)
  y.fill <- y.fill |> group_by(plot,symbol,taxon,type,stratum,stratum.label, stratum.min, stratum.max,Bottom,Top,dbh.Low, dbh.High, cover.ps) |> summarise(Cover = max(Cover), BA = max(BA), frq.strat=max(frq.strat))
  #get quantiles in consideration of zeros for absences
  y.fill <- y.fill |> group_by(taxon, symbol,type, stratum, stratum.label, stratum.min, stratum.max, Bottom, Top, dbh.Low, dbh.High, cover.ps) |>
    summarise(cover.Low = round(quantile(Cover, lowerQ),1),
              cover.mean = mean(Cover),
              cover.High = round(quantile(Cover, upperQ),1),
              BA.Low = round(quantile(BA, lowerQ),1),
              BA.mean = mean(BA),
              BA.High = round(quantile(BA, upperQ),1),
              frq.strat = round(mean(frq.strat),3))
  y.fill <- y.fill |> group_by(symbol, taxon, type) |> mutate(taxon.cover = cover.agg(cover.mean), over.cover = cover.agg(ifelse(Top > 5,cover.mean,0)))
  y.fill <- y.fill |> group_by(type) |> mutate(type.top = max(Top))
  y.fill <- left_join(y.fill, f) |> mutate(BA.pp = round(BA.mean/frq.plot,1),
                                           cover.pp = round(cover.mean/frq.plot,1),
                                           BA.mean = round(BA.mean,1),
                                           cover.mean = round(cover.mean,1),
                                           frq.plot = round(frq.plot,3))
  
  y.fill <- subset(y.fill, select=c(taxon, symbol, type, stratum, stratum.label, stratum.min, stratum.max, cover.Low, cover.High, Bottom,Top, dbh.Low, dbh.High, BA.Low, BA.High, cover.mean, cover.pp, cover.ps, BA.mean, BA.pp, frq.plot, frq.strat, taxon.cover,over.cover, type.top)) |> arrange(-type.top, type, -over.cover, -taxon.cover, -Top)
  
  
  return(y.fill)
}




breaks = c(0.5,2,5,12,24)
veg.summ <- multi.releve.summary(veg1, breaks)
veg.summ0 <- summary.ESIS0(veg1, breaks)
veg.con <- consolidated.summary(veg.summ, breaks)
write.csv(veg.con, 'veg.con.csv', row.names = FALSE, na = "")




consolidated.summary <- function(veg.summ, breaks = c(0.5,5,12)){
  
  if(length(breaks) < 4){
    cnames = c(paste0('0-',breaks[1],'m'),paste0(breaks[1],'-',breaks[2],'m'),
               paste0(breaks[2],'-',breaks[3],'m'),paste0(breaks[3],'+ m'),'ht max (m)')
    veg.con <- veg.summ |> group_by(habitsort, overstory, covertotal, taxon, type, freq) |> summarise(c1 = sum(ifelse(stratum %in% 1, cover.r,0)),
                                                                                                      c2 = sum(ifelse(stratum %in% 2, cover.r,0)),
                                                                                                      c3 = sum(ifelse(stratum %in% 3, cover.r,0)),
                                                                                                      c4 = sum(ifelse(stratum %in% 4, cover.r,0)),
                                                                                                      ht.max = max(ht.max)) |>
      arrange(habitsort, -overstory, -covertotal, taxon) |> mutate(#habitsort=NULL, overstory=NULL, covertotal=NULL,
        c1=ifelse(c1<=0,NA,c1),
        c2=ifelse(c2<=0,NA,c2),
        c3=ifelse(c3<=0,NA,c3),
        c4=ifelse(c4<=0,NA,c4))
    colnames(veg.con)[colnames(veg.con) %in% c('c1','c2','c3','c4','ht.max')] <- cnames
    
    
  }else{
    cnames = c(paste0('0-',breaks[1],'m'),paste0(breaks[1],'-',breaks[2],'m'),
               paste0(breaks[2],'-',breaks[3],'m'),
               paste0(breaks[3],'-',breaks[4],'m'),paste0(breaks[4],'+ m'),'ht max (m)')
    
    veg.con <- veg.summ |> group_by(habitsort, overstory, covertotal, taxon, type, freq) |> summarise(
      c1 = sum(ifelse(stratum %in% 1, cover.r,0)),
      c2 = sum(ifelse(stratum %in% 2, cover.r,0)),
      c3 = sum(ifelse(stratum %in% 3, cover.r,0)),
      c4 = sum(ifelse(stratum %in% 4, cover.r,0)),
      c5 = sum(ifelse(stratum %in% 5, cover.r,0)),
      ht.max = max(ht.max)) |>
      arrange(habitsort, -overstory, -covertotal, taxon) |> mutate(#habitsort=NULL, overstory=NULL, covertotal=NULL,
        c1=ifelse(c1<=0,NA,c1),
        c2=ifelse(c2<=0,NA,c2),
        c3=ifelse(c3<=0,NA,c3),
        c4=ifelse(c4<=0,NA,c4),
        c5=ifelse(c5<=0,NA,c5))
    colnames(veg.con)[colnames(veg.con) %in% c('c1','c2','c3','c4','c5','ht.max')] <- cnames}
  return(veg.con)}





veg.summ <- summary.ESIS0(veg1, breaks)





flat.summary <- function(veg.summ){
  breaks <- sort(unique(veg.summ$stratum.min)) 
  breaks <- breaks[breaks>0]
  nbks <- length(breaks)+1
  brks <- c(0,breaks,1000)
  veg.con <- veg.summ |> group_by(taxon,symbol,type,frq.plot, taxon.cover,over.cover,type.top) |> summarise(
    ht.max = max(Top))
  for(i in 1:nbks){#i=6
    y <- veg.summ |> group_by(taxon,symbol,type,frq.plot, taxon.cover,over.cover,type.top) |> 
      summarise(c0 = sum(ifelse(stratum %in% i, cover.mean,0))) |> mutate(c0 = ifelse(c0 > 0, c0, NA_real_))
    veg.con <- veg.con |> cbind(c0=y$c0)
    colnames(veg.con)[colnames(veg.con) %in% 'c0'] <- paste0('c',i)
    cname <- paste0('c',i)
    if(i==1){cnames=cname}else{cnames=c(cnames,cname)}
    
    altname <- case_when(i == 1  ~ paste0('0-',brks[i],' m'),
                       i < nbks ~ paste0(brks[i],'-',brks[i+1],' m'),
                       TRUE ~ paste0(brks[i],'+ m'))
    if(i==1){altnames=altname}else{altnames=c(altnames,altname)}
  }
  veg.con <- veg.con |> 
      arrange(-type.top, type, -over.cover, -taxon.cover, -ht.max, taxon) |> subset(select = c('taxon','type','frq.plot',cnames,'ht.max','over.cover','taxon.cover'))
    colnames(veg.con)[colnames(veg.con) %in% c(cnames,'frq.plot')] <- c('frq',altnames)
  return(veg.con)}


breaks = c(0.5,2,5,12,24)


veg.summ <- summary.ESIS0(veg1, breaks)
veg.con <- flat.summary(veg.summ)








