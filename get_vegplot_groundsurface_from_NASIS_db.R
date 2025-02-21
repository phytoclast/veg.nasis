

#' Get vegetation plot transect ground surface cover data from local NASIS database
#'
#' @param SS fetch data from the currently loaded selected set in NASIS or from the entire local database (default: TRUE)
#' @param dsn Optional: path to local SQLite database containing NASIS table structure; default: NULL
#' @param si Convert legacy units to SI (above ground measurements); default: TRUE
#'
#' @return Data frame containing a summary data for line point intercept ground surface cover hits by cover type.
#' @export
#'
#' @examples vsurf <- get_vegplot_groundsurface_from_NASIS_db(SS=F)
#' @examples allplots <- vsurf[,c("siteiid","siteobsiid","vegplotid","vegplotname","vegtransectid","totalpointssampledcount","transectlength")] |> unique()
#' @examples alltypes <- vsurf[,c("groundsurfcovtype","groundcoverptcount","groundcoverptpct")] |> mutate(groundcoverptcount=0, groundcoverptpct=0) |> unique()
#' @examples vmerge <- merge(allplots, alltypes)
#' @examples vsurf <- vsurf |> select(colnames(vmerge)) |> rbind(vmerge)
#' @examples vsurf <- vsurf |> group_by(vegplotid, vegplotname, groundsurfcovtype) |> summarise(cover = sum(groundcoverptcount/totalpointssampledcount*100))
#' @examples vsurf <- vsurf |> mutate(treament=substr(vegplotname,1,4))
#' @examples vsurfsum <- vsurf |> group_by(treament, groundsurfcovtype) |> summarise(cover05 = quantile(cover, 0.05),
#' @examples                                                                         covermean = mean(cover),
#' @examples                                                                        cover95 = quantile(cover, 0.95))
get_vegplot_groundsurface_from_NASIS_db <- function(SS = TRUE, dsn = NULL, si = TRUE) {
  q <- "SELECT siteiid, siteobsiid, vegplotid, vegplotname, vegtransectid, vt.totalpointssampledcount, vt.transectlength, groundsurfcovtype, groundcoverptcount, groundcoverptpct, quadratsize, quadratshape, groundcoverquadpctave
              FROM site_View_1 AS s
              INNER JOIN siteobs_View_1 AS so ON so.siteiidref=s.siteiid
              INNER JOIN vegplot_View_1 AS v ON v.siteobsiidref=so.siteobsiid
              INNER JOIN vegtransect_View_1 AS vt
                     ON vt.vegplotiidref = v.vegplotiid
              LEFT JOIN transectgroundsurfcover_View_1 AS vtps
                     ON vtps.vegtransectiidref = vt.vegtransectiid"
              # LEFT JOIN groundsurfcovdetails_View_1 AS vtpsd
              #        ON vtpsd.transectgrsurfcoviidref = vtps.transectgroundsurfcoveriid

  if (!SS) {
    q <- gsub("_View_1", "", q)
  }

  res <- dbQueryNASIS(NASIS(dsn = dsn), q)
  if(si){res$transectlength <- round(res$transectlength*0.3048,0)}
  uncode(res)
}






#' Get vegetation plot transect plant foliar cover data from local NASIS database
#'
#' @param SS fetch data from the currently loaded selected set in NASIS or from the entire local database (default: TRUE)
#' @param dsn Optional: path to local SQLite database containing NASIS table structure; default: NULL
#' @param si Convert legacy units to SI (above ground measurements); default: TRUE
#'
#' @return Data frame containing a summary data for line point intercept foliar cover hits by plant taxon.
#'
#' @export
#'
#' @examples vfoliar <- get_vegplot_transplantfoliar_from_NASIS_db(SS=F)
get_vegplot_transplantfoliar_from_NASIS_db <- function(SS = TRUE, dsn = NULL, si=TRUE) {
  q <- "SELECT siteiid, siteobsiid, vegplotid, vegplotname, vt.vegtransectid, vt.totalpointssampledcount,
  vt.transectlength, plantsym, plantsciname, plantnatvernm, speciesfoliarcovhitcount, speciesfoliarcovpctlineint,plantheightcllowerlimit,plantheightclupperlimit
              FROM site_View_1 AS s
              INNER JOIN siteobs_View_1 AS so ON so.siteiidref=s.siteiid
              INNER JOIN vegplot_View_1 AS v ON v.siteobsiidref=so.siteobsiid
              INNER JOIN vegtransect_View_1 AS vt
                     ON vt.vegplotiidref = v.vegplotiid
              LEFT JOIN vegtransectplantsummary_View_1 AS vtps
                     ON vtps.vegtransectiidref = vt.vegtransectiid
              LEFT JOIN plant ON plant.plantiid = vtps.plantiidref"
  if (!SS) {
    q <- gsub("_View_1", "", q)
  }

  res <- dbQueryNASIS(NASIS(dsn = dsn), q)
  if(si){
    res$transectlength <- round(res$transectlength*0.3048,0)
    res$plantheightcllowerlimit <- round(res$plantheightcllowerlimit*0.3048,1)
    res$plantheightclupperlimit <- round(res$plantheightclupperlimit*0.3048,1)
  }
  uncode(res)
}

#' Get vegetation plot snag and basal area from local NASIS database
#'
#' @param SS fetch data from the currently loaded selected set in NASIS or from the entire local database (default: TRUE)
#' @param dsn Optional: path to local SQLite database containing NASIS table structure; default: NULL
#' @param si Convert legacy units to SI (above ground measurements); default: TRUE
#'
#' @return
#' @export Data frame containing a plot summary of snag density and tree basal area,
#'
#' @examples vplot <- get_vegplot_BA_snags_from_NASIS_db(SS=T, si=FALSE) |> mutate(treatment=substr(vegplotname,1,4))
#' @examples vplot <- vplot |> group_by(treatment) |> summarise(snagsoft_low = min(treesnagdensitysoft), snagsoft_high = max(treesnagdensitysoft),       snaghard_low = min(treesnagdensityhard), snaghard_high = max(treesnagdensityhard))

get_vegplot_BA_snags_from_NASIS_db <- function(SS = TRUE, dsn = NULL, si=TRUE) {
  q <- "SELECT siteiid, siteobsiid, vegplotid, vegplotname, overstorycancontotalpct, basalareaplottotal, treesnagdensityhard, treesnagdensitysoft
              FROM site_View_1 AS s
              INNER JOIN siteobs_View_1 AS so ON so.siteiidref=s.siteiid
              INNER JOIN vegplot_View_1 AS v ON v.siteobsiidref=so.siteobsiid
              "
  if (!SS) {
    q <- gsub("_View_1", "", q)
  }

  res <- dbQueryNASIS(NASIS(dsn = dsn), q)
  if(si){
    res$basalareaplottotal <- round(res$basalareaplottotal*0.229568,1)
    res$treesnagdensityhard <- round(res$treesnagdensityhard/0.4046856,0)
    res$treesnagdensitysoft <- round(res$treesnagdensitysoft/0.4046856,0)
  }
  uncode(res)
}
