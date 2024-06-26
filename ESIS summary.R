
library(vegnasis)
veg.raw <- vegnasis::nasis.veg
veg <- clean.veg(veg.raw) |> fill.type.df() |> fill.hts.df()
veg <- veg |> mutate(type= fill.type.ESIS(taxon))
veg <- veg |> mutate(symbol = fill.usda.symbols(taxon, symbol))
veg <- subset(veg,  grepl('2022MI165',plot))

veg.esis <- summary.ESIS(veg, lowerQ = 0.5, upperQ = 0.75)



x=veg
summary.ESIS <-  function(x, breaks=c(0.5,5,15), lowerQ=0.25, upperQ=0.75){
  x=x |> mutate(dbh.min= ifelse(is.na(dbh.min), dbh.max,dbh.min))
  #frequency of whole plot
f <- x |> subset(cover>0, select=c(plot,taxon)) |> unique()  |> mutate(freq=1)
nplots = length(unique(f$plot))
f <- f |> group_by(taxon) |> summarise(plot.frq=round(sum(freq)/nplots*100,1))

  y <- NULL
  nbks <- length(breaks)+1
  brks <- c(0,breaks,1000)
  #extract means by stratum and plot
  for(i in 1:(nbks)){#i = 8
    y0 <- x %>% subset(ht.max < brks[i+1] & ht.max >= brks[i])

    if(nrow(y0)>0){
      y0 <- y0 %>% mutate(stratum=i, stratum.label = paste0(brks[i], "-", ifelse(i==nbks, "+",brks[i+1])), ht.min= ht.min, ht.max = ht.max)
      y1 <- y0 %>% group_by(plot, symbol, taxon, type, stratum, stratum.label) %>%
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
  y = y |> group_by(symbol,taxon,type,stratum, stratum.label) |>
    mutate(Bottom=round(weighted.mean(ht.min, Cover+0.001, na.rm=TRUE),1),
           Top=round(weighted.mean(ht.max, Cover+0.001, na.rm=TRUE),1),
           dbh.Low =  weighted.mean(dbh.min, Cover+0.001, na.rm=TRUE),
           dbh.High =  weighted.mean(dbh.max, Cover+0.001, na.rm=TRUE),
           dbh.Low =  ifelse(is.nan(dbh.Low), NA, round(dbh.Low,0)),
           dbh.High =  ifelse(is.nan(dbh.High), NA, round(dbh.High,0)),
           cover.present = round(mean(Cover),1))

  y = y |> group_by(plot) |> mutate(totalBA = sum(BA, na.rm = TRUE), overCover = ifelse(Top > 5, Cover, NA), grossCover = sum(overCover, na.rm = TRUE), BA = round(totalBA*Cover/(grossCover+0.000001),1))
#get frequency in stratum
  y = y |> group_by(plot,symbol,taxon,type) |> mutate(frq = ifelse(sum(Cover)>0,1,0))
 #insert zeros for missing species found in other plots
  y.plot <- unique(subset(y, select=c("plot")))
  y.mid <- unique(subset(y, select=c("symbol","taxon","type","stratum","stratum.label","Bottom","Top","dbh.Low","dbh.High","cover.present")))
  y.fill <- merge(y.plot, y.mid) |> mutate(Cover = 0, BA = 0, frq=0)
  y.fill <- rbind(y, y.fill)
  y.fill <- y.fill |> group_by(plot,symbol,taxon,type,stratum,stratum.label,Bottom,Top,dbh.Low, dbh.High, cover.present) |> summarise(Cover = max(Cover), BA = max(BA), frq=max(frq))
#get quantiles in consideration of zeros for absences
  y.fill <- y.fill |> group_by(taxon, symbol,type, Bottom,Top, dbh.Low, dbh.High, cover.present) |>
    summarise(cover.Low = round(quantile(Cover, lowerQ),1),
              cover.mean = round(mean(Cover),1),
              cover.High = round(quantile(Cover, upperQ),1),
              BA.Low = round(quantile(BA, lowerQ),1),
              BA.High = round(quantile(BA, upperQ),1),
              frq = round(mean(frq)*100,1))
  y.fill <- y.fill |> group_by(symbol, taxon, type) |> mutate(taxon.cover = cover.agg(cover.mean))
  y.fill <- y.fill |> group_by(type) |> mutate(type.top = max(Top))
  y.fill <- left_join(y.fill, f)

  y.fill <- subset(y.fill, select=c(taxon, symbol, type, cover.Low, cover.High, Bottom,Top, dbh.Low, dbh.High, BA.Low, BA.High, cover.mean, cover.present, frq, plot.frq, taxon.cover, type.top)) |> arrange(-type.top, type, -taxon.cover, -Top)


  return(y.fill)
}




veg.str <- summary.strata(veg,  breaks=c(0.5,5,15)) |> structure.fill.zero() |> subset(type %in% c('tree', 'shrub/vine', 'grass/grasslike',  'forb'))
