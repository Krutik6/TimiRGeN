#' @title getP
#' @description Internal function that retreived p value from a linear
#' regression model.
#' @param modelobject linear regression model that has been generated by the
#' linearRegr function.
#' @return
#' @noRd
#' @importFrom stats pf
getP <- function (modelobject) {

  if (!is(modelobject, "lm"))  stop("Not an object of class 'lm' ")

  f <- summary(modelobject)$fstatistic

  p <- pf(f[1],f[2],f[3],lower.tail=FALSE)

  attributes(p) <- NULL

  return(p)

}
