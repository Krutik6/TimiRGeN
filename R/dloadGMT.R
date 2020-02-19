#' @title dloadGMT
#' @description Provides mouse or human wikipathway information.
#' @param MAE MultiAssayExperiment to store downlaoded GMT data in.
#' @param speciesInitials Either Hs or Mm, respectively to retreive Homo
#'sapiens or Mus musculus data.
#' @return 3 dataframes. 1) genes-wikipathwayIDs 2) wikipathwayIDs-wikipathway
#'fullnames 3) combination of 1 and 2. All of which can be sotred in an MAE.
#' @export
#' @importFrom clusterProfiler read.gmt
#' @usage dloadGMT(MAE, speciesInitials = "")
#' @examples
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' dloadGMT(MAE, speciesInitial = "Mm") -> MAE
dloadGMT <- function(MAE, speciesInitials){
if (missing(speciesInitials)) stop('speciesInitials should be either Mm
for mouse data or Hs for human data.')
if (speciesInitials == "Mm") {
download.file(paste("http://data.wikipathways.org/20200110/gmt/",
"wikipathways-20200110-gmt-Mus_musculus.gmt", sep = ""),
"mus.gmt")
gmt <- read.gmt("mus.gmt")
file.remove("mus.gmt")
} else if (speciesInitials == "Hs") {
download.file(paste("http://data.wikipathways.org/20200110/gmt/",
"wikipathways-20200110-gmt-Homo_sapiens.gmt", sep = ""),
"hom.gmt")
gmt <- read.gmt("hom.gmt")
file.remove("hom.gmt")
}
ont <- wpid <- gene <- name <- NULL
pathways <- gmt %>% tidyr::separate(ont, c("name","version","wpid",
"org"), "%")
MAE@ExperimentList$path_gene <- as.data.frame(pathways%>%dplyr::select(wpid, 
gene))
MAE@ExperimentList$path_name <- as.data.frame(pathways%>%dplyr::select(wpid, 
name))
MAE@ExperimentList$path_data <- as.data.frame(pathways%>%dplyr::select(wpid,
gene,
name))
return(MAE)
}
