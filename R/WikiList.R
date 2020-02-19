#' @title WikiList
#' @description Provides a list of wikipathways and each gene found to be
#' within that wikipathway as the value.
#' @details install rWikiPathways.
#' @param stringSpecies species wikipathways relate to e.g. 'Homo sapiens'
#' @param stringSymbol Type of gene code e.g. 'H'
#' @return List of wikipathways and associated genes as a list of strings.
#' @export
#' @import rWikiPathways
#' @usage WikiList(stringSpecies = '', stringSymbol = '')
WikiList <- function(stringSpecies, stringSymbol){
if (missing(stringSpecies)) stop('Input a spieces name e.g. Homo
sapiens or Mus musculus.');
if (missing(stringSymbol)) stop('Input a symbol type e.g. H, L or
En.');
pathlist <- rWikiPathways::listPathways(stringSpecies)
lapply(pathlist, `[[`, 1) -> wpslist
lapply(wpslist, function(x) rWikiPathways::getXrefList(x,
stringSymbol)) -> listolists
names(listolists) <- wpslist
return(listolists)
}

