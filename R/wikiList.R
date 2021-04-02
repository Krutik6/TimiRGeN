#' @title wikiList
#' @description Provides a list of wikipathways specific for a species, and each
#' gene found to be within each pathway. wikiList will download and process a
#' large amount of data from the wikipathways website so this may take some
#' time to complete.
#' @param MAE MultiAssayExperiment which will store the results from
#' wikiList. It is recommended to use the same MAE which stores output from the
#' dloadGmt function.
#' @param stringSpecies Full species name to decide which wikipathways to
#' download. 'Homo sapiens' to download human pathways or 'Mus musculus'
#' to download mouse pathways.
#' @param stringSymbol Type of gene ID to retrieve e.g. 'En' for ensemble gene
#' IDs or 'L' for entrezgene IDs.
#' @return List of wikipathways and associated genes saved as as strings. Output
#' will be stored as metadata in the input MAE.
#' @export
#' @importFrom rWikiPathways listPathways getXrefList
#' @usage wikiList(MAE, stringSpecies = '', stringSymbol = '')
wikiList <- function(MAE, stringSpecies, stringSymbol){

    if (missing(MAE)) stop('MAE is missing. Add MAE. Results of wikiList will be stored in this MAE.')

    if (missing(stringSpecies)) stop('stringSpecies is missing. Add a species name e.g. "Homo sapiens" or "Mus musculus".')

    if (missing(stringSymbol)) stop('stringSymbol is missing. Add a symbol type e.g. "En" or "L".')

    metadata <- `metadata<-` <- NULL

    # get all wikipathways for a species
    pathlist <- rWikiPathways::listPathways(stringSpecies)

    # Convert column of IDs into list
    y <- as.list(pathlist$id)

    names(y) <- pathlist$id

    # Find gene names for each pathways
    wikilist <- lapply(y, function(x) rWikiPathways::getXrefList(x,
                                                            stringSymbol))

    metadata(MAE)[["wikilist"]] <- wikilist

return(MAE)
}

