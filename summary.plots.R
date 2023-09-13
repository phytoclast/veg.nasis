setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
veg.spp <- read.delim('data/Observed_Species.txt')
veg.site <- read.delim('data/Sites.txt')
veg <- clean.veg.log(veg.site,veg.spp)

veg <- veg |> fill.type.df() |> fill.hts.df()
groups <- c('2023OH051003', '2023OH051001','2023OH051002', '2023523.001', '2023523.002')
veg1 <- subset(veg, label %in% groups)



x = veg1
breaks=c(0.5,5,15)

multi.releve.summary <-  function(x, breaks=c(0.5,5,15), woodytypes = c('tree','shrub/vine', 'epiphyte')){
  
  y <- NULL
  nbks <- length(breaks)+1
  brks <- c(0,breaks,1000)
  
  
  for(i in 1:(nbks)){#i = 1
    y0 <- x %>% subset((ht.max < brks[i+1] & ht.max >= brks[i] & type %in% woodytypes)|(i==1 & !type %in% woodytypes))
    
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
  y <- y |> left_join(freq0) |> mutate(cover.p = round(cover/freq,1), cover.r = round(cover,1))
  y <- subset(y, select=c(taxon, type, stratum, stratum.label, freq, cover.r, cover.p, bottom, top, ht.min, ht.max))
  y <- y |> mutate(habitsort = case_when(type %in% 'tree' ~ 1,
                                         type %in% 'shrub/vine' ~ 2,
                                         type %in% 'grass/grasslike' ~ 3,
                                         type %in% 'forb' ~ 4,
                                         TRUE ~ 5))
  
  y <- y |> group_by(taxon, type) |> mutate(overstory = cover.agg(ifelse(top > 5,1,0)*cover.r), covertotal = cover.agg(cover.r)) |> arrange(habitsort, -overstory, -covertotal, taxon, -top)
                                         
  
  return(y)
}



veg.summ <- multi.releve.summary(veg1, c(0.5,5,12))

