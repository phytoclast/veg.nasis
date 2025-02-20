get_vegplot_groundsurface_from_NASIS_db <- function(SS = TRUE, dsn = NULL) {
  q <- "SELECT vegplotid, vegplotname, vt.totalpointssampledcount, vt.transectlength,  vegtransectiidref, groundsurfcovtype, groundcoverptcount, groundcoverptpct, quadratsize, quadratshape, groundcoverquadpctave, transectgroundsurfcoveriid
              FROM site_View_1 AS s
              INNER JOIN siteobs_View_1 AS so ON so.siteiidref=s.siteiid
              INNER JOIN vegplot_View_1 AS v ON v.siteobsiidref=so.siteobsiid
              INNER JOIN vegtransect_View_1 AS vt
                     ON vt.vegplotiidref = v.vegplotiid
              LEFT JOIN transectgroundsurfcover_View_1 AS vtps
                     ON vtps.vegtransectiidref = vt.vegtransectiid
"
  if (!SS) {
    q <- gsub("_View_1", "", q)
  }

  res <- dbQueryNASIS(NASIS(dsn = dsn), q)
  res$transectlength_m <- round(res$transectlength*0.3048,0)
  uncode(res)
}
transurf <- get_vegplot_groundsurface_from_NASIS_db(SS=F)
