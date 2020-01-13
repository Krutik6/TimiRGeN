#' @title dloadmiRDB
#' @description Downloads most recent version (6.0) of predicted targets from
#'targetscan database.
#' @return A dataframe of predicted miR-mRNA interactions in different
#' species.
#' @export
#' @usage dloadmiRDB()
dloadmiRDB <- function(){
download.file("http://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
'MIRDB.gz')
read.table(gzfile("MIRDB.gz"), header = FALSE, fill = TRUE) -> miRDB
file.remove("MIRDB.gz")
print("Downloaded Predicted Targets from miRDB version 6.0")
print("Check if a newer version is available? http://mirdb.org/." )
print("If so download that and run it through miRDB_data function.")
return(miRDB)
}
