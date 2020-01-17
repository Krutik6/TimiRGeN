#' @title downloadGMT
#' @description Provides mouse or human wikipathway information.
#' @param speciesInitials Either Hs or Mm, respectively to retreive Homo
#'sapiens or Mus musculus data.
#' @return 3 dataframes. 1) genes-wikipathwayIDs 2) wikipathwayIDs-wikipathway
#'fullnames 3) combination of 1 and 2.
#' @export
#' @importFrom clusterProfiler read.gmt
#' @usage downloadGMT(speciesInitials = "")
#' @examples
#'downloadGMT(speciesInitial = "Mm")
#'file.remove("mus.gmt")
downloadGMT <- function(speciesInitials){
if (missing(speciesInitials)) stop('speciesInitials should be either Mm
for mouse data or Hs for human data.')
if (speciesInitials == "Mm") {
download.file(paste("http://data.wikipathways.org/20190610/gmt/",
"wikipathways-20190610-gmt-Mus_musculus.gmt", sep = ""),
"mus.gmt")
read.gmt("mus.gmt") -> gmt
file.remove("mus.gmt")
} else if (speciesInitials == "Hs") {
download.file(paste("http://data.wikipathways.org/20190610/gmt/",
"wikipathways-20190610-gmt-Homo_sapiens.gmt", sep = ""),
"hom.gmt")
ont <- wpid <- gene <- name <- NULL
read.gmt("hom.gmt") -> gmt
file.remove("hom.gmt")
}
pathways <- gmt %>% tidyr::separate(ont, c("name","version","wpid",
"org"), "%")
path_gene <- pathways %>% dplyr::select(wpid, gene)
path_name <- pathways %>% dplyr::select(wpid, name)
path_data <- pathways %>% dplyr::select(wpid, gene, name)
path_list <- list(path_gene = path_gene,
path_name = path_name,
path_data = path_data)
return(list2env(path_list, .GlobalEnv))
}
