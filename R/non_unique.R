#' @title non_unique
#' @description Internal function for getIDs_miR_mousetohuman. Helps data
#' wrangling for -3p and -5p genes.
#' @param Col Column of dataframe to modify.
#' @param sep Seperator to insert.
#' @param suffix Suffix to add behind seperator, default is 1, and increases
#' in numerical order for each repeated string in the column.
#' @return Modified column of a dataframe. Each non-unique value will have
#' numbers representing it's duplicate number behing it.
#' @usage non_unique(Col, sep, suffix)
#' @export
#' @examples
#'data.frame(name = c("mmu-miR-101a-3p", "mmu-miR-101a-5p", "mmu-miR-101c",
#'"mmu-miR-106a-5p", "mmu-miR-106b-3p"),
#'Hs_n = c("hsa-miR-101-1", "hsa-miR-101-1", "hsa-miR-101c",
#'"hsa-miR-106a", "hsa-miR-106b")) -> sample
#'non_unique(Col = sample$Hs_n, sep = '-', suffix = 'p') -> sample$Hs_n
non_unique <- function(Col, sep, suffix){
ave(as.character(Col), Col,
FUN=function(x) if (length(x)>1
) paste0(x[1], sep, seq_along(x), suffix) else x[1]) -> Col
return(Col)
}
