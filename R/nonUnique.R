#' @title nonUnique
#' @description Helps to create the adjusted miR ID dataframes during the
#' getIDsMir functions. This is an internal function.
#' @param Col Column of dataframe to modify. sep and suffix will be added behind
#' repeated strings within this column..
#' @param sep Separator to insert. Will separate suffix from current string in
#' selected column.
#' @param suffix Suffix to be added behind the sep parameter, default is 1,
#' and it increases in numerical order for each repeated string in the column.
#' @return Modifies a column of a dataframe. Each non-unique value will have
#' a unique number added behind it.
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
