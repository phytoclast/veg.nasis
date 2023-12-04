setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
library(sf)
library(terra)
library(soilDB)
library(aqp)
library(stringi)
plants <- readRDS('data/plants/plants20231129.RDS')
#saveRDS(plants, 'D:/scripts/veg.nasis/data/plants/plants20231129.RDS')

kew <- read.table('data/plants/wcvp__2_/wcvp_names.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
geo <- read.table('data/plants/wcvp__2_/wcvp_distribution.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

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
    enc = stringi::stri_enc_detect(x[i])
    enc <- stringi::enc[[1]]$Encoding[1]
    if(!enc %in% 'UTF-8'){
      if(enc %in% 'windows-1252'){
        x[i] <- stringi::stri_conv(x[i], from = 'windows-1252', to='UTF-8')
      }else if(enc %in% 'ISO-8859-1'){
        x[i] <- stringi::stri_conv(x[i], from = 'ISO-8859-1', to='UTF-8')
      }else{
        x[i] <- stringi::stri_conv(x[i], from = 'latin1', to='UTF-8')
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

m.ac <- read.csv('data/plants/m.ac.csv', encoding = 'UTF-8')
#m.ac <- m.ac |> mutate(en1 = whatEncoding(syn), en2 = whatEncoding(syn.auth))

m.ac <- m.ac |> mutate(acc = cleanEncoding(acc),
                       acc.auth = cleanEncoding(acc.auth),
                       syn = cleanEncoding(syn),
                       syn.auth = cleanEncoding(syn.auth))
m.ac <- m.ac |> mutate(acc = extractTaxon(acc),
                       syn = extractTaxon(syn))


#create new three way synonymy
#bonap
syn <- unique(data.frame(taxon = m.ac$syn, author = m.ac$syn.auth, bonap = m.ac$acc))
#kew
wfo <- unique(data.frame(id = kew$plant_name_id,ac = kew$accepted_plant_name_id, taxon = kew$taxon_name, author=kew$taxon_authors, status = kew$taxon_status))
wfo <- subset(wfo, status %in% c("Synonym","Accepted","Orthographic"))
wfo <- wfo |> mutate(taxon = extractTaxon(taxon))
wfo.acc <- wfo[,c("id","taxon")]; colnames(wfo.acc)<-c("ac","kew")
wfo <- wfo |> left_join(wfo.acc)
wfo <- wfo[,c("taxon", "kew")]
syn <- syn |> left_join(wfo, multiple = 'first')
colnames(wfo) <- c("bonap", "kew2" )
syn <- syn |> left_join(wfo, multiple = 'first')
syn <- syn |> mutate(kew = ifelse(is.na(kew), kew2,kew), kew2 = NULL)
#usda
plants <- subset(plants, !auct)
plants <- plants |> mutate(Synonym.Symbol = ifelse(is.na(Synonym.Symbol)|Synonym.Symbol %in% "", Symbol,Synonym.Symbol))
plants.acc <- plants[,c("Symbol","taxon")]
colnames(plants.acc) <- c("Symbol","usda")
plants <- plants |> left_join(plants.acc, multiple = 'first')
syn <- syn |> left_join(plants[,c('taxon','usda')], multiple = 'first')
plants$usda2 <- plants$usda
syn <- syn |> left_join(plants[,c('taxon','usda2')], by=join_by(bonap==taxon), multiple = 'first')
syn <- syn |> mutate(usda = ifelse(is.na(usda), usda2,usda), usda2 = NULL)
#recycle bonap synonyms 1. kew
syn.bonapkew <- subset(syn, !is.na(kew), select = c(bonap, kew))
colnames(syn.bonapkew) <- c("bonap","kew2")
syn <- left_join(syn, syn.bonapkew, multiple = 'first')
syn <- syn |> mutate(kew = ifelse(is.na(kew), kew2,kew), kew2 = NULL)
#recycle bonap synonyms 2. usda
syn.bonapusda <- subset(syn, !is.na(usda), select = c(bonap, usda))
colnames(syn.bonapusda) <- c("bonap","usda2")
syn <- left_join(syn, syn.bonapusda, multiple = 'first')
syn <- syn |> mutate(usda = ifelse(is.na(usda), usda2,usda), usda2 = NULL)
write.csv(syn, 'data/plants/syns.csv', row.names = F)


# syn <- syn |>  mutate(infrataxon0 = extractTaxon(taxon, 'infrataxon'),
#                       epithet1 = extractTaxon(bonap, 'epithet'),
#                       infrataxon1 = extractTaxon(bonap, 'infrataxon'), 
#                       binomial2 = extractTaxon(kew, 'binomial'),
#                       epithet2 = extractTaxon(kew, 'epithet'),
#                       infrataxon2 = extractTaxon(kew, 'infrataxon'),
#                       binomial3 = extractTaxon(usda, 'binomial'),
#                       epithet3 = extractTaxon(usda, 'epithet'),
#                       infrataxon3 = extractTaxon(usda, 'infrataxon'))
# syn1 <- subset(syn, infrataxon0 == '' & infrataxon1 == '' & epithet3 == infrataxon3 & !is.na(usda) & grepl(' ',usda))


syns <- vegnasis::syns



bonap.missing <- subset(bonap, !bonap %in% syn$taxon)


orthvar <-  subset(plants, grepl('orth. var.', author))
orthvar <-  subset(plants, Symbol  %in% orthvar$Symbol)


variants <- as.data.frame(rbind(
  c('j','j'),
  c('i','j'),
  c('io','o'),
  c('o','io'),
  c('i','y'),
  c('y','i'),
  c('ck','k'),
  c('k','ck'),
  c('iae','ica'),
  c('ica','iae'),
  c('pseudo','pseudo-'),
  c('pseudo-','pseudo'),
  c('nn','n'),
  c('n','nn'),
  c('a$','us'),
  c('us$','a'),
  c('ensis$','ense'),
  c('ense$','ensis'),
  c('[^i^e]i$','ii'),
  c('ii$','i'),
  c('ae','i'),
  c('i','ae'),
  c('ii','ae'),
  c('a$','um'),
  c('um$','a'),
  c('e$','is'),
  c('is$','e'),
  c('ped','paed'),
  c('paed','ped'),
  c('ces','caes'),
  c('caes','ces'),
  c('aea$','ea'),
  c('[^a]ea$','aea'),
  c('r$','ra'),
  c('ra$','r'),
  c('r$','rus'),
  c('rus$','r'),
  c('r$','rum'),
  c('rum$','r')
));colnames(variants) <- c('from','to')
