# load the package
# devtools::install_github("jinyizju/V.PhyloMaker")
# devtools::install_github("jinyizju/V.PhyloMaker2")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("V.PhyloMaker2")
# input example species list
c1 <- c("Araucaria araucana", "Sequoiadendron giganteum", "Ephedra antisyphilitica", "Cycas revoluta", "Ginkgo biloba", "Callicarpa americana","Cornus florida","Rhododendron maximum","Sassafras albidum","Asimina triloba", "Magnolia grandiflora", "Liatris spicata", "Eupatorium altissimum", "Ageratina altissima", "Eutrochium maculatum","Carya ovata", "Quercus palustris", "Acer rubrum", "Platycarya strobilacea", "Tilia americana", "Zea mays","Schizachyrium scoparium","Poa pratensis", "Sabal palmetto","Trillium grandiflorum","Liriodendron tulipifera", "Amborella trichopoda","Pinus palustris","Osmunda regalis","Lycopodium clavatum")
c2 <- c("Araucaria", "Sequoiadendron", "Ephedra", "Cycas", "Ginkgo","Callicarpa","Cornus","Rhododendron", "Sassafras","Asimina", "Magnolia", "Liatris","Eupatorium", "Ageratina", "Eutrochium","Carya", "Quercus", "Acer", "Platycarya", "Tilia", "Zea","Schizachyrium","Poa", "Sabal","Trillium","Liriodendron","Amborella","Pinus","Osmunda","Lycopodium")
c3 <- c("Araucariaceae", "Cupressaceae", "Ephedraceae", "Cycadaceae", "Ginkgo","Lamiaceae","Cornaceae","Ericaceae", "Lauraceae","Annonaceae", "Magnoliaceae", "Asteraceae","Asteraceae", "Asteraceae", "Asteraceae","Juglandaceae", "Fagaceae", "Sapindaceae", "Juglandaceae", "Malvaceae", "Poaceae","Poaceae","Poaceae", "Arecaceae","Melanthiaceae", "Magnoliaceae","Amborellaceae", "Pinaceae","Osmundaceae", "Lycopodiaceae")
example <- data.frame(species = c1, genus = c2, family = c3)
taxa.TPL = tips.info.TPL
taxa.WP = tips.info.WP
taxa.LCVP = tips.info.LCVP
library(dplyr)
apg <- read.csv('data/plants/apg.csv')
apgselect <- subset(apg, superclass %in% c("Acrogymnosperms"))
taxa.select <- subset(taxa.TPL, family %in% apgselect$family)
taxa.select <- taxa.select |> slice_sample(n=3, by=family, replace = TRUE) |> unique()
example <- unique(data.frame(species = taxa.select$species, genus = taxa.select$genus, family = taxa.select$family))

### run the function
result <- phylo.maker(example)

### plot the phylogenies with node ages displayed.

w <- 800
h <- nrow(example)*12+80
u <- 12
png(filename='phyogeny.png',width = w, height = h, units = "px", pointsize = u)

par(mar = c(2,0,1,13))
plot.phylo(result$scenario.3, cex = 1.5, main = "Phylogeny with dates")
nodelabels(round(branching.times(result$scenario.3), 1), cex = 1)

dev.off()

