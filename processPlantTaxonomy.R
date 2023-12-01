setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
library(sf)
library(terra)
library(soilDB)
library(aqp)
library(stringi)
plants <- read.csv('D:/scripts/veg.nasis/data/plants/USDA Plants Database20231129.txt')
saveRDS(plants, 'D:/scripts/veg.nasis/data/plants/plants20231129.RDS')

kew <- read.table('D:/scripts/veg.nasis/data/plants/wcvp__2_/wcvp_names.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 
geo <- read.table('D:/scripts/veg.nasis/data/plants/wcvp__2_/wcvp_distribution.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")  

geo.select <- subset(geo, continent %in% 'NORTHERN AMERICA' | area %in% c('Puerto Rico','Hawaii'))
kew.select <- subset(kew, plant_name_id %in% geo.select$plant_name_id)
kew.select2 <- subset(kew, accepted_plant_name_id %in% kew.select$plant_name_id)



plants <- plants |> mutate(taxon = extractTaxon(Scientific.Name.with.Author, 'taxon'),
                           author = extractTaxon(Scientific.Name.with.Author, 'author'))
plants <- plants |> mutate(auct = ifelse(grepl('auct',author) |
                                           grepl('illeg.',author) |
                                           grepl(' non',author), T,F))

cleanEncoding <- function(x){
  n = length(x)
  for(i in 1:n){#i=1
    enc = stri_enc_detect(x[i])
    enc <- enc[[1]]$Encoding[1]
    if(!enc %in% 'UTF-8'){
      if(enc %in% 'windows-1252'){
        x[i] <- stri_conv(x[i], from = 'windows-1252', to='UTF-8')
      }else if(enc %in% 'ISO-8859-1'){
        x[i] <- stri_conv(x[i], from = 'ISO-8859-1', to='UTF-8')
      }else{
        x[i] <- stri_conv(x[i], from = 'latin1', to='UTF-8')
      }
    }
  }  
  return(x)}


whatEncoding <- function(x){
  n = length(x)
  for(i in 1:n){#i=1
    enc <- stri_enc_detect(x[i])
    enc <- enc[[1]]$Encoding[1]
    x[i] <- enc}
  return(x)}

m.ac <- read.csv('D:/scripts/vegnasis/data_raw/m.ac.csv', encoding = 'UTF-8')
#m.ac <- m.ac |> mutate(en1 = whatEncoding(syn), en2 = whatEncoding(syn.auth))

m.ac <- m.ac |> mutate(acc = cleanEncoding(acc),
                       acc.auth = cleanEncoding(acc.auth),
                       syn = cleanEncoding(syn), 
                       syn.auth = cleanEncoding(syn.auth))
write.csv(m.ac, 'D:/scripts/vegnasis/data_raw/m.ac.csv', row.names = F)

stri_conv(this, from = enc, to='utf-8')



x = m.ac.latin1[57957,'syn.auth']
stri_conv(x, from = 'windows-1254', to='UTF-8')
stri_conv(x, from = 'windows-1250', to='UTF-8')

enc <- stri_enc_detect(m.ac.latin1[57957,'syn.auth'])
enc <- enc[[1]]$Encoding[1]
enc = t(as.data.frame(enc))
enc[[1]]$Encoding[1]

m.ac.utf8 <- read.csv('D:/scripts/vegnasis/data_raw/m.ac.csv', encoding = 'UTF-8')
m.ac.latin1 <- read.csv('D:/scripts/vegnasis/data_raw/m.ac.csv', encoding = 'latin1')
m.ac <- m.ac |> mutate(acc = extractTaxon(acc), syn = extractTaxon(syn))

bonap <- unique(data.frame(bonap = m.ac$syn, bonap.auth = m.ac$syn.auth))
bonap <- bonap

syn <- unique(data.frame(taxon = kew.select2$taxon_name, author=kew.select2$taxon_authors))

syn <- syn |> mutate(taxon = extractTaxon(taxon))

syn 
