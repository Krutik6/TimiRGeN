#' @title dloadMirdb
#' @description Downloads most recent version (6.0) of predicted targets from
#'targetscan database.
#' @param MAE MultiAssayExperiment object.
#' @return A dataframe of predicted miR-mRNA interactions in different
#' species.
#' @export
#' @usage dloadMirdb(MAE)
#' @examples
#' \dontrun{
#' library(MultiAssayExperiment)
#' MAE <- MultiAssayExperiment()
#' MAE <-dloadMirdb(MAE)
#' }
dloadMirdb <- function(MAE){
    if (missing(MAE)) stop('Add MultiAssayExperiment');
    download.file("http://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
                  'MIRDB.gz')

    miRDB <- read.table(gzfile("MIRDB.gz"), header = FALSE, fill = TRUE)
    file.remove("MIRDB.gz")

    MAE2 <- suppressMessages(MultiAssayExperiment(list('miRDB' = miRDB)))
    MAE <- c(MAE, MAE2)

    print("Downloaded Predicted Targets from miRDB version 6.0")
    print("Check if a newer version is available? http://mirdb.org/." )
    print("If so download that and run it through miRDB_data function.")

return(MAE)
}
