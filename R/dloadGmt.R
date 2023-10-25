#' @title dloadGmt
#' @description Downloads the most up-to-date versions of the mouse or human
#' wikipathway information databases. Output will be stored as three distinct
#' dataframes within the input MAE 1) path_gene, 2) path_names, 3) path_data.
#' @param MAE MultiAssayExperiment to store downloaded GMT data in. It might
#' be useful to start a new MAE for dloadGmt using MultiAssayExperiment(). This
#' is so the MAE objects used in this analysis do not get too large.
#' @param species Full species names e.g. retrieve "Homo sapiens" or "Mus musculus" data.
#' @return 3 dataframes. 1) path_gene, 2) path_names, 3) path_data.
#' All of which will be stored as assays in the input MAE.
#' @export
#' @importFrom clusterProfiler read.gmt
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @importFrom stringr %>%
#' @usage dloadGmt(MAE, species = "")
#' @examples
#' MAE <- MultiAssayExperiment()
#'
#' MAE <- dloadGmt(MAE, species = "Homo sapiens")
dloadGmt <- function(MAE, species){

  if (missing(MAE)) stop('MAE is missing. Add MultiAssayExperiment to store files downloaded by dloadGmt. The user may wish to create a new MAE object using the MultiAssayExperiment() function.')

  if (missing(species)) stop('species is missing. E.g. "Rattus norvegicus" for rat data and "Homo sapiens" for human data.')

  ont <- wpid <- gene <- name <- NULL

  # Download data
  t <- try(x <- rWikiPathways::downloadPathwayArchive(organism=species,
                                                      format = "gmt"))

  #back-up download if current gmt files not yet ready
  if("try-error" %in% class(t)){x <- rWikiPathways::downloadPathwayArchive(
                                                      date = "20230310",
                                                      organism=species,
                                                      format = "gmt")}
  # load file
  gmt <- clusterProfiler::read.gmt(x)

  # delete file
  file.remove(x)

  gmt <- as.data.frame(gmt)

  colnames(gmt) <- c("ont", "gene")

  # separate by % into organised data frame
  pathways <- gmt %>% tidyr::separate(ont, c(
                                  "name","version","wpid","org"),
                                    "%")

  # create genes, names and data files from gmt file
  path_gene <- as.data.frame(pathways%>%dplyr::select(wpid, gene))

  path_name <- as.data.frame(pathways%>%dplyr::select(wpid, name))

  path_data <- as.data.frame(pathways%>%dplyr::select(wpid,gene,name))

  MAE2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list("path_gene" = path_gene,
                                                     "path_name" = path_name,
                                                     "path_data" = path_data))
  ))

  MAE <- suppressWarnings(c(MAE, MAE2))

  return(MAE)
}
