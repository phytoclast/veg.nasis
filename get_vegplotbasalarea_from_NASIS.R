library(soilDB)


get_vegplot_speciesbasalarea_from_NASIS <- function(SS = TRUE, dsn = NULL) {
  q <- "SELECT siteiid, siteobsiid, vegplotiid, vegplotid, vegplotname, obsdate, 
  primarydatacollector, plantiidref AS plantiid, 
  plotspeciebasalareaiid, basalareatreescountediid
            plantsym, plantsciname, plantnatvernm,
            basalareafactor, speciesnumbertreesin, speciesbasalarea,
            treenumber, treeheight, treediameterbreastheight
FROM site_View_1 AS s
  INNER JOIN siteobs_View_1 AS so ON so.siteiidref = s.siteiid
  LEFT JOIN vegplot_View_1 AS v ON v.siteobsiidref = so.siteobsiid
  LEFT JOIN plotspeciesbasalarea_View_1 AS vb ON vb.vegplotiidref = v.vegplotiid
    LEFT JOIN basalareatreescounted_View_1 AS ba  ON ba.plotspeciebasalareaiidref = vb.plotspeciebasalareaiid
    INNER JOIN plant ON plant.plantiid = vb.plantiidref"
  
  if (!SS) {
    q <- gsub("_View_1", "", q)
  }
  
  uncode(dbQueryNASIS(NASIS(dsn = dsn), q), dsn = dsn)
}


veg.basal.raw <- get_vegplot_speciesbasalarea_from_NASIS(SS=F)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

saveRDS(veg.basal.raw, 'data/veg.basal.raw.RDS')