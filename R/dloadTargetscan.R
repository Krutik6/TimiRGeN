#' @title dloadTargetscan
#' @description Downloads most recent version (7.2) of predicted targets from
#' targetscan database.
#' @param MAE MultiAssayExperiment object.
#' @return A dataframe of predicted miR-mRNA interactions in different
#' species.
#' @export
#' @usage dloadTargetscan(MAE)
#' @examples
#' \dontrun{
#'
#' MAE <- MultiAssayExperiment()
#'
#' MAE <-dloadTargetscan(MAE)
#' }
dloadTargetscan <- function(MAE){

    if (missing(MAE)) stop('Add MultiAssayExperiment.')

    # download file
    download.file(paste("http://www.targetscan.org/vert_72/vert_72_data_download/",
                        "Predicted_Targets_Context_Scores.default_predictions.txt.zip",
                         sep = ''), 'Targetscans.zip')

    # unzip and then remove zipped file
    unzip(zipfile = 'Targetscans.zip')

    file.remove("Targetscans.zip")

    # read unzipped file
    TargetScans <- read.table('Predicted_Targets_Context_Scores.default_predictions.txt',
                              fill = TRUE, header = TRUE)

    TargetScans <- TargetScans[,seq_len(5)]

    # remove unzipped file
    file.remove("Predicted_Targets_Context_Scores.default_predictions.txt")

    MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                                  'TargetScans' = TargetScans)))
    MAE <- c(MAE, MAE2)

    # Message on download version
    print("Downloaded Predicted Targets Context Scores default version 7.2")

    print("Check if a newer version is available?")

    print("http://www.targetscan.org/vert_72/")

    print("If so download and run it through Targetscans_data function.")

return(MAE)
}
