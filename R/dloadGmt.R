#' @title dloadGmt
#' @description Provides mouse or human wikipathway information.
#' @param MAE MultiAssayExperiment to store downlaoded GMT data in.
#' @param speciesInitials Either Hs or Mm, respectively to retreive Homo
#'sapiens or Mus musculus data.
#' @return 3 dataframes. 1) genes-wikipathwayIDs 2) wikipathwayIDs-wikipathway
#'fullnames 3) combination of 1 and 2. All of which can be sotred in an MAE.
#' @export
#' @importFrom clusterProfiler read.gmt
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @usage dloadGmt(MAE, speciesInitials = "")
#' @examples
#' miR <- mm_miR
#' mRNA <- mm_mRNA
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
dloadGmt <- function(MAE, speciesInitials){

    if (missing(MAE)) stop('Add MultiAssayExperiment');
    if (missing(speciesInitials)) stop('speciesInitials should be either Mm
                                       for mouse data or Hs for human data.');

    ont <- wpid <- gene <- name <- NULL

    # If mouse is selected?
    if (speciesInitials == "Mm") {
        download.file(paste("http://data.wikipathways.org/20200110/gmt/",
                            "wikipathways-20200110-gmt-Mus_musculus.gmt",
                            sep = ""),
                            "mus.gmt")
        gmt <- read.gmt("mus.gmt")
        file.remove("mus.gmt")

    # If human is selected?
    } else if (speciesInitials == "Hs") {
        download.file(paste("http://data.wikipathways.org/20200110/gmt/",
                            "wikipathways-20200110-gmt-Homo_sapiens.gmt",
                            sep = ""),
                            "hom.gmt")
        gmt <- read.gmt("hom.gmt")
        file.remove("hom.gmt")
    }

    gmt <- as.data.frame(gmt)
    colnames(gmt) <- c("ont", "gene")
    # separate by % into data frame
    pathways <- gmt %>% separate(ont, c("name","version","wpid","org"), "%")

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
