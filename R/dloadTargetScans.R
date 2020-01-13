#' @title dloadTargetScans
#' @description Downloads most recent version (7.2) of predicted targets from
#' targetscan database.
#' @return A dataframe of predicted miR-mRNA interactions in different
#' species.
#' @export
#' @usage dloadTargetScans()
#' @examples dloadTargetScans() -> TargetScans
dloadTargetScans <- function(){
download.file(paste("http://www.targetscan.org/vert_72/vert_72_data_download/",
"Predicted_Targets_Context_Scores.default_predictions.txt.zip",
sep = ''), 'Targetscans.zip')
unzip(zipfile = 'Targetscans.zip')
print("Downloaded Predicted Targets Context Scores default version 7.2")
print("Check if a newer version is available?")
print("http://www.targetscan.org/vert_72/")
print("If so download that and run it through Targetscans_data function.")
file.remove("Targetscans.zip")
TargetScans <- read.table(
'Predicted_Targets_Context_Scores.default_predictions.txt',
fill = TRUE, header = TRUE)
TargetScans[,seq_len(5)] -> TargetScans
return(TargetScans)
}
