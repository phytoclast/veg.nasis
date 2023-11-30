setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
library(sf)
library(terra)
library(soilDB)
library(aqp)

plants <- read.csv('D:/scripts/veg.nasis/data/plants/USDA Plants Database20231129.txt')
saveRDS(plants, 'D:/scripts/veg.nasis/data/plants/plants20231129.RDS')

kew <- read.table('D:/scripts/veg.nasis/data/plants/wcvp__2_/wcvp_names.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 
geo <- read.table('D:/scripts/veg.nasis/data/plants/wcvp__2_/wcvp_distribution.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")  

geo.select <- subset(geo, continent %in% 'NORTHERN AMERICA' | area %in% c('Puerto Rico','Hawaii'))
kew.select <- subset(kew, plant_name_id %in% geo.select$plant_name_id)
kew.select2 <- subset(kew, accepted_plant_name_id %in% kew.select$plant_name_id)

m.ac <- read.csv('D:/scripts/vegnasis/data_raw/m.ac.csv')

plants <- read.csv('D:/scripts/veg.nasis/data/plants/USDA Plants Database20231129.txt')
extractTaxon <- function(rawnames, report = 'taxon'){
x = data.frame(rawname = rawnames)
x <- x |> mutate(rawname = str_replace_all(rawname,'^x\\s|^X\\s','×'),
                           rawname = str_replace_all(rawname,'\\sx\\s|\\sX\\s',' ×'),
                           genus = trimws(str_split_fixed(rawname,'[[:space:]]',2)[,1]),
                           first = trimws(str_split_fixed(rawname,'[[:space:]]',2)[,2]),
                           second = ifelse(grepl('^[a-z]|^×',first) & !grepl('^ex\\w',first),
                                           trimws(str_split_fixed(first,'[[:space:]]',2)[,1]),""),
                           xxx = ifelse(grepl('^[a-z]|^×',first) & !grepl('^ex\\w',first),
                                        trimws(str_split_fixed(first,'[[:space:]]',2)[,2]),first),
                           xxx = str_replace_all(xxx, "\\sssp.\\s|^ssp.\\s|\\ssubsp. |^subsp.\\s", '_ssp._' ),
                           xxx = str_replace_all(xxx, "\\svar.\\s|^var.\\s", '_var._' ),
                           xxx = ifelse(grepl('\\sf.\\s[a-z]',xxx)|grepl('^f.\\s[a-z]',xxx) & !grepl('_f._ex\\s',xxx), str_replace_all(xxx, "\\sf.\\s|^f.\\s", '_f._' ),xxx),                         xxx1 = trimws(str_split_fixed(xxx,'_',3)[,1]),
                           xxx2 = trimws(str_split_fixed(xxx,'_',3)[,2]),
                           xxx3 = trimws(str_split_fixed(xxx,'_',4)[,3]),
                           xxx4 = trimws(str_split_fixed(xxx,'_',5)[,4]),
                           xxx5 = trimws(str_split_fixed(xxx,'_',6)[,5]),
                           epithet = second,
                           binomialauthor = xxx1,
                           infrarank1 = xxx2,
                           infraname1 = trimws(str_split_fixed(xxx3,'[[:space:]]',2)[,1]),
                           infraauthor1 = trimws(str_split_fixed(xxx3,'[[:space:]]',2)[,2]),
                           infrarank2 = xxx4,
                           infraname2 = trimws(str_split_fixed(xxx5,'[[:space:]]',2)[,1]),
                           infraauthor2 = trimws(str_split_fixed(xxx5,'[[:space:]]',2)[,2]),
                           # first=NULL,
                           # second=NULL,
                           # xxx=NULL,
                           # xxx1=NULL,
                           # xxx2=NULL,
                           # xxx3=NULL,
                           # xxx4=NULL,
                           # xxx5=NULL,
                           binomial = trimws(paste(genus, epithet)),
                           taxon = trimws(paste(binomial, infrarank1, infraname1, infrarank2, infraname2)),
                           author = trimws(ifelse(!infraauthor2 %in% "", infraauthor2, ifelse(!infraauthor1 %in% "", infraauthor1, binomialauthor)))
)
if(report == 'author'){return(x$author)}else if(report=='binomial'){return(x$binomial)}else if(report == 'genus'){return(x$genus)}else{return(x$taxon)}
}



rawnames = c('Arnica angustifolia Vahl subsp. tomentosa (Macoun) G.W. Douglas & G. Ruyle-Douglas', 'Abies ×shastensis (Lemmon) Lemmon', 'Abies X shastensis (Lemmon) Lemmon','Abies balsamea ssp. lasiocarpa', 'Abies balsamea (L.) Mill. var. fallax (Engelm.) B. Boivin','Agastache pallidiflora (A. Heller) Rydb. ssp. pallidiflora var. greenei (Briq.) R.W. Sanders')
extractTaxon(rawnames,'taxon')