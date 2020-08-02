#' @title wikiList
#' @description Provides a list of wikipathways and each gene found to be
#' within that wikipathway. This will be stored as metadata of a MAE. wikiList
#' will download a large amount of data from the wikipathways website so this
#' may take some time.
#' @param MAE MultiAssayExperiment object which the results of wikiList is
#' to be added. It is recommended to use the same MAE which stores data from
#' dloadGmt.
#' @param stringSpecies Species name to decide which wikipathways to download.
#' 'Homo sapiens' to download human pathways or 'Mus musculus' to download
#' mouse pathways.
#' @param stringSymbol Type of gene code e.g. 'H' for human gene names
#' , 'En' for ensemble gene IDs , 'L' for entrezgene IDs.
#' @return List of wikipathways and associated genes as a list of strings.
#' @export
#' @importFrom rWikiPathways listPathways getXrefList
#' @usage wikiList(MAE, stringSpecies = '', stringSymbol = '')
#' @examples
#' \donttest{
#' MAE <- MultiAssayExperiment()
#'
#' MAE <- wikiList(MAE, stringSpecies = 'Homo sapiens', stringSymbol = 'En')
#' }
wikiList <- function(MAE, stringSpecies, stringSymbol){

    if (missing(MAE)) stop('Add MAE object. Results of wikiList will be
                           stored in this MAE.')

    if (missing(stringSpecies)) stop('Input a spieces name e.g. "Homo
                                      sapiens" or "Mus musculus".')

    if (missing(stringSymbol)) stop('Input a symbol type e.g. "H", "L" or "En."
                                    ')

    metadata <- `metadata<-` <- NULL

    # get all wikipathways for a species
    pathlist <- rWikiPathways::listPathways(stringSpecies)

    wpslist <- lapply(pathlist, `[[`, 1)
    # Find gene names for each pathways
    wikilist <- lapply(wpslist, function(x) rWikiPathways::getXrefList(x,
                                                            stringSymbol))
    names(wikilist) <- wpslist

    metadata(MAE)[["wikilist"]] <- wikilist

return(MAE)
}

