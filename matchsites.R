library(vegnasis)
set.site <- function(plot, lat, lon, date = NA_integer_, year = NA_integer_, month = NA_integer_, day = NA_integer_){
  x <- data.frame(plot = plot,
                  lat = lat,
                  lon = lon,
                  date = date,
                  year = year,
                  month = month,
                  day = day
  )
  x <- x |> mutate(date = as.Date(ifelse(is.na(date),as.character(as.Date(ISOdate(year = year,
                                                                                  month = month,
                                                                                  day = day))),as.character(date)))
                   ,
                   # year = ifelse(is.na(year), as.integer(format(date, format="%Y")),as.integer(year)),
                   # month = ifelse(is.na(month), as.integer(format(date, format="%m")),as.integer(month)),
                   # day = ifelse(is.na(day), as.integer(format(date, format="%d")),as.integer(day)),
                   year = NULL,
                   month = NULL,
                   day = NULL)
  return(x)
}

match.sites <- function(x,y,maxdist = 30, maxdays = NA){

  xplots <- x$plot
  x$link <- NA_character_

  for(i in 1:length(xplots)){#i=500
    yf <- y
    if(!is.na(maxdays)){
      yf <- subset(yf, abs(date - x[i,]$date)<=maxdays)
    }

    if(nrow(yf) > 0){
      yf$dist <- (((x[i,]$lat - yf$lat)/360*40041.47*1000)^2 +
                    ((x[i,]$lon - yf$lon)/360*40041.47*1000*cos(x[i,]$lat/2/360*2*3.141592))^2)^0.5

      mindist <- min(yf$dist, na.rm = TRUE)
      y0 <- yf %>% subset(dist %in% mindist)
      x[i,] <- x[i,] |> mutate(link = ifelse(mindist <= maxdist,  y0$plot[1], link))
    }
  }
  return(x)
}



x <-vegnasis::obs
y <- soilDB::get_site_data_from_NASIS_db(SS=F)




x <- set.site(plot=x$Observation_Label,
              lat = x$Latitude,
              lon = x$Longitude,
              year = x$Year,
              month = x$Mon,
              day= x$Day)
y <- set.site(plot=y$site_id,
              lat = y$y,
              lon = y$x,
              date = y$obs_date)

x <- x |> subset(year == 2022 & month == 06 & day %in% c(22,23))

clean.veg()
z = match.sites(x=x,y=y, maxdist = 10000, maxdays = 200)
