#' @title dloadMirtarbase
#' @description Downloads most recent version (7.0) of functional targets from
#'miRTarBase database.
#' @param MAE MultiAssayExperiment object.
#' @return Dataframe of functionally tested miR-mRNA interactions in different
#' species.
#' @export
#' @usage dloadMirtarbase(MAE)
#' @importFrom readxl read_excel
#' @examples
#' \dontrun{
#'
#' MAE <- MultiAssayExperiment()
#'
#' MAE <-dloadMirdb(MAE)
#'
#' }
dloadMirtarbase <- function(MAE){

  if (missing(MAE)) stop('Add MultiAssayExperiment.')

  # download miRTarBase data
  download.file(
  'http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/miRTarBase_MTI.xlsx',
  'miRTarBase.xlsx')

  # read using readxl
  miRTarBase <- readxl::read_excel('miRTarBase.xlsx')

  # write it as a CSV file so it can be imported back in the intended structure
  write.csv(miRTarBase, 'miRTarBase.csv')

  miRTarBase <- read.csv('miRTarBase.csv', row.names = 1)

  # Remove multiple files
  file.remove("miRTarBase.xlsx")

  file.remove("miRTarBase.csv")

  # add to MAE object
  MAE2 <- suppressMessages(MultiAssayExperiment(list('miRTarBase' = miRTarBase)
                                                                              ))

  MAE <- c(MAE, MAE2)

  # Information message
  print("Downloaded published miRNA target interaction data version 7.0")

  print("Check if a newer version is available?")

  print("http://mirtarbase.mbc.nctu.edu.tw/php/index.php")

  print("If so download that and run it through miRTarBase_data function.")

return(MAE)
}
