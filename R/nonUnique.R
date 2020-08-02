#' @title nonUnique
#' @description Internal function for getIDs functions. Helps data
#' wrangling for -3p and -5p genes.
#' @param Col Column of dataframe to modify.
#' @param sep Seperator to insert. Will separate suffix from current string in
#' column.
#' @param suffix Suffix to add behind separator, default is 1, and increases
#' in numerical order for each repeated string in the column.
#' @return Modified column of a dataframe. Each non-unique value will have
#' numbers representing it's duplicate number behind it.
#' @importFrom stats ave
#' @usage nonUnique(Col, sep, suffix)
#' @noRd
nonUnique <- function(Col, sep, suffix){

    # If string is seen more than once add a suffix to it
    Col <- stats::ave(as.character(Col), Col,
                FUN=function(x) if (length(x)>1) paste0(x[1],
                                                        sep,
                                                        seq_along(x),
                                                        suffix) else x[1])
return(Col)
}
