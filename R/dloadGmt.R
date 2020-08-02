#' @title dloadGmt
#' @description Provides mouse or human wikipathway information. This data is
#' downloaded from http://data.wikipathways.org/20200110/gmt/. This will be
#' stored as three distinct data frames within a MAE object
#' 1) path_gene 2) path_names 3) path_data. It is recommended that the
#' user creates a new MAE object to avoid MAE objects getting too large.
#' @param MAE MultiAssayExperiment to store downloaded GMT data in. It might
#' be useful to start a new MAE for dloadGmt using MultiAssayExperiment(). This
#' is so the MAE objects used in this analysis do not get too large.
#' @param speciesInitials Either "Hs" or "Mm", respectively to retrieve Homo
#'sapiens or Mus musculus data.
#' @return 3 dataframes. 1) path_gene 2) path_names 3) path_data.
#' All of which will be stored in the input MAE, in the assay section.
#' @export
#' @importFrom clusterProfiler read.gmt
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @importFrom stringr %>%
#' @usage dloadGmt(MAE, speciesInitials = "")
#' @examples
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- MultiAssayExperiment()
#'
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
dloadGmt <- function(MAE, speciesInitials){

    if (missing(MAE)) stop('Add MultiAssayExperiment to store files downloaded
                            by dloadGmt. The user may wish to create a new
                            MAE object using the MultiAssayExperiment()
                            function.')

    if (missing(speciesInitials)) stop('speciesInitials should be either "Mm"
                                       for mouse data or "Hs" for human data.')

    ont <- wpid <- gene <- name <- NULL

    # If mouse is selected, download mouse gmt data to file
    if (speciesInitials == "Mm") {
        download.file(paste("http://data.wikipathways.org/20200110/gmt/",
                            "wikipathways-20200110-gmt-Mus_musculus.gmt",
                            sep = ""),
                            "mus.gmt")

        # load file
        gmt <- clusterProfiler::read.gmt("mus.gmt")

        # delete file
        file.remove("mus.gmt")

    # If human is selected, download human gmt data from file
    } else if (speciesInitials == "Hs") {
        download.file(paste("http://data.wikipathways.org/20200110/gmt/",
                            "wikipathways-20200110-gmt-Homo_sapiens.gmt",
                            sep = ""),
                            "hom.gmt")

        # load file
        gmt <- read.gmt("hom.gmt")

        # delete file
        file.remove("hom.gmt")
    }

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

    MAE2 <- suppressMessages(MultiAssayExperiment(list("path_gene" = path_gene,
                                                       "path_name" = path_name,
                                                       "path_data" = path_data)
                                                  ))

    MAE <- c(MAE, MAE2)

    return(MAE)
}
