
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(vegnasis)
familylink <- read.csv('data/plants/familylink.csv')
m.ac <- read.csv('data/plants/m.ac.csv')

gensyn <- m.ac |> subset(select = c(ac.genus, sy.genus)) |> unique()
gensyn <- gensyn |> left_join(familylink, by=c('ac.genus'='genus'))

synfam <- gensyn |> subset(select = c(sy.genus, family)) |> unique()
synfam <- synfam |> group_by(sy.genus) |> mutate(ct = length(family))
synfam <- subset(synfam, ct == 1 & !sy.genus %in% familylink$genus & !is.na(family))
synfam$genus = synfam$sy.genus
synfam$ac <- FALSE
familylink$ac <- TRUE
familylink.comb <- familylink |> rbind(synfam[,colnames(familylink)])

write.csv(familylink.comb, 'C:/workspace2/vegnasis/data_raw/familylink.csv', row.names = FALSE)
familylink <- read.csv('C:/workspace2/vegnasis/data_raw/familylink.csv')
usethis::use_data(familylink, overwrite = T)

apg <- read.csv('C:/workspace2/veg.nasis/data/plants/apg.csv')
write.csv(apg, 'C:/workspace2/vegnasis/data_raw/apg.csv', row.names = FALSE)
apg <- read.csv('C:/workspace2/vegnasis/data_raw/apg.csv')
usethis::use_data(apg, overwrite = T)
