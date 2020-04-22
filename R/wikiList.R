#' @title wikiList
#' @description Provides a list of wikipathways and each gene found to be
#' within that wikipathway as the value.
#' @details install rWikiPathways.
#' @param MAE MultiAssayExperiment object.
#' @param stringSpecies species wikipathways relate to e.g. 'Homo sapiens'
#' @param stringSymbol Type of gene code e.g. 'H'
#' @return List of wikipathways and associated genes as a list of strings.
#' @export
#' @import rWikiPathways
#' @usage wikiList(MAE, stringSpecies = '', stringSymbol = '')
#' @examples
#' \donttest{
#' library(MultiAssayExperiment)
#' MAE <- MultiAssayExperiment()
#' MAE <- wikiList(MAE, stringSpecies = 'Homo sapiens', stringSymbol = 'En')
#' }
wikiList <- function(MAE, stringSpecies, stringSymbol){

    if (missing(MAE)) stop('Add MAE object');
    if (missing(stringSpecies)) stop('Input a spieces name e.g. Homo
    sapiens or Mus musculus.');
    if (missing(stringSymbol)) stop('Input a symbol type e.g. H, L or
    En.');

    metadata <- `metadata<-` <- NULL

    pathlist <- rWikiPathways::listPathways(stringSpecies)
    wpslist <- lapply(pathlist, `[[`, 1)
    wikilist <- lapply(wpslist, function(x) rWikiPathways::getXrefList(x,
                                                            stringSymbol))
    names(wikilist) <- wpslist
    metadata(MAE)[["wikilist"]] <- wikilist
return(MAE)
}

