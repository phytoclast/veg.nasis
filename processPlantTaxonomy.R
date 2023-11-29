plants <- read.csv('D:/scripts/veg.nasis/data/plants/USDA Plants Database20231129.txt')
saveRDS(plants, 'D:/scripts/veg.nasis/data/plants/plants20231129.RDS')

kew <- read.table('D:/scripts/veg.nasis/data/plants/wcvp__2_/wcvp_names.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 
geo <- read.table('D:/scripts/veg.nasis/data/plants/wcvp__2_/wcvp_distribution.csv', sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")  

geo.select <- subset(geo, continent %in% 'NORTHERN AMERICA' | area %in% c('Puerto Rico','Hawaii'))
kew.select <- subset(kew, plant_name_id %in% geo.select$plant_name_id)
kew.select2 <- subset(kew, accepted_plant_name_id %in% kew.select$plant_name_id)