

ytrans = 'identity'; yratio=1; units = 'm'; skycolor = "#D9F2FF80"; fadecolor = "#D9F2FF"; gridalpha=0.3; groundcolor="#808066"; xlim=c(0,50); ylim=c(-1, 30+5); xticks=5; yticks=5; xslope=15; yslope=25; xperiod=10; xamplitude=0


library(vegnasis)
veg.raw <- vegnasis::nasis.veg
veg <- clean.veg(veg.raw) |> fill.type.df() |> fill.hts.df()
veg$habit <- get.habit.code(veg$taxon)
veg.s <- subset(veg,  grepl('2022MI165021.P',plot))
taxon <- c('Acer rubrum', 'Pinus resinosa')
crfill <- c(NA,"#80991A")
stfill <- c('gray',"#B36666")
crshape <- c(NA,'conifer2')
override <- data.frame(taxon=taxon,stfill=stfill,crfill=crfill,crshape=crshape)
veg.s <- veg.s |> left_join(override)

timeA=Sys.time()
veg_profile_plot0(plants, xslope=50, yslope=25)
Sys.time()-timeA
timeA=Sys.time()
veg_profile_plot(plants)
Sys.time()-timeA
plants <- grow_plants(veg.s, pwidth=50)

#test looping ggplot, test color layering; test any other way to group ggplot by yp
xnmax <- max(plants$xn, na.rm =TRUE)
xnmin <- min(plants$xn, na.rm =TRUE)
ypmax <- max(plants$yp, na.rm =TRUE)
ypmin <- min(plants$yp, na.rm =TRUE)
ypwid <- ypmax-ypmin
plants <- plants |> arrange(yp,stumpid, objid, ptord) |> mutate(zn = zn+(xp*xslope/100)+((yp-ypmin)*yslope/100)+
                                                                  xamplitude+xamplitude*sin(xp/xperiod*3.141592*2)) #implement slope
zmax <- max(plants$zn, na.rm =TRUE)

# colormixer('magenta', fadecolor, round(1-1/(1+((18.19-ypmin)/20)),2))

iters = 20
ypmin <- min(plants$yp)-0.01
ypmax <- max(plants$yp)+0.01
ypinc <- (ypmax - ypmin)/iters

# 1-(floor(((ypmax - yp)/(ypmax - ypmin))*iters+0.99)/iters))
# yp=ypmax - c(1:20)*ypinc
# 1-(floor(((ypmax - yp)/(ypmax - ypmin))*5+0.99)/5)
plants <- plants |> mutate(fill=colormixer(fill, fadecolor, round(1-1/(1+((yp-ypmin)/20)),2)), 
                           color=colormixer(color, fadecolor, round(1-1/(1+((yp-ypmin)/20)),2)))
pcolor <- c(plants$color) |> unique() |> sort()
pfill <- c(plants$fill) |> unique()|> sort()
ucf = case_when(units %in% c('feet', 'ft') ~ 0.3048,
                units %in% c('inches', 'in') ~ 0.3048/12,
                units %in% c('cm') ~ 0.01,
                TRUE ~ 1)
units = ifelse(ucf == 1, 'm',units)

yunits = paste0('height (', units,')')
xunits = paste0('ground distance (', units,')')
ybreaks = seq(floor(ylim[1]/ucf/yticks)*yticks-yticks,
              floor(ylim[2]/ucf/yticks)*yticks+yticks,
              yticks)*ucf
xbreaks = seq(floor(xlim[1]/ucf/xticks)*xticks-xticks,floor(xlim[2]/ucf/xticks)*xticks+xticks,xticks)*ucf
yminor = seq(floor(ylim[1]/ucf-yticks),floor(ylim[2]/ucf+yticks),yticks/5)*ucf
xminor = seq(floor(xlim[1]/ucf-xticks),floor(xlim[2]/ucf+xticks),xticks/5)*ucf
ylabels = ybreaks/ucf
xlabels =  xbreaks/ucf




gp <- ggplot()
for(i in 1:iters){#i=1 i=1:iters
  # fade = 1-(floor(((i-1)/iters)*5)/5+0.2)
  ypcut0 = ypmax - i*ypinc
  ypcut1 = ypmax - (i-1)*ypinc
plants0 <- plants |> subset(yp > ypcut0 & yp <= ypcut1)
stems0 <- plants0 |> subset(obj %in% 'stem')
  crowns0 <- plants0 |> subset(obj %in% c('crown','herb'))
  gp = gp+
    geom_polygon(data=stems0, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)+
    geom_polygon(data=crowns0, aes(x=xn,y=zn,group=objid, fill=fill, color=color), alpha=1, linewidth=0.01)
}

gp <- gp+
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
  scale_y_continuous(name = yunits, trans = ytrans, labels = ylabels, breaks = ybreaks, minor_breaks = yminor, limits = c(-10,zmax+5))+#
  scale_x_continuous(name = xunits ,breaks = xbreaks, labels = xlabels, minor_breaks = xminor, limits = c(xnmin-5,xnmax+5))#

gp






p=0.46

colormixer <- function(colorname, mixcolor, p){
  ccc <- col2rgb(colorname)
  ccc <- data.frame(r = ccc[1,],   g = ccc[2,],   b = ccc[3,])
  mmm <- col2rgb(mixcolor)
  new <- ccc |> mutate(r = r*(1-p)+mmm[1,1]*p,
                       g = g*(1-p)+mmm[2,1]*p,
                       b = b*(1-p)+mmm[3,1]*p)
  new <- rgb(new$r,new$g,new$b, maxColorValue = 255)
  return(new)
}





