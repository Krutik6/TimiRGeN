#' @title dloadMirtarbase
#' @description Downloads most recent version (7.0) of functional targets from
#'miRTarBase database.
#' @param MAE MultiAssayExperiment object.
#' @param species mmu or hsa.
#' @return Dataframe of functionally tested miR-mRNA interactions in different
#' species.
#' @export
#' @usage dloadMirtarbase(MAE, species)
#' @examples
#' library(MultiAssayExperiment)
#' MAE <- MultiAssayExperiment
#' MAE <- dloadMirtarbase(MAE, species = "hsa")
dloadMirtarbase <- function(MAE, species){

    if (missing(MAE)) stop('Add MultiAssayExperiment');
    if (missing(species)) stop('Select mmu to hsa');

    mmu_miRTarBase <- hsa_miRTarBase <- NULL

    if (species == "mmu") {
        miRTarBase <- mmu_miRTarBase
        MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                                      'miRTarBase' = miRTarBase
                                                      )))
        MAE <- c(MAE, MAE2)
        return(MAE)
    } else if (species == "hsa") {
        miRTarBase <- hsa_miRTarBase
        MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                                      'miRTarBase' = miRTarBase
                                                      )))
        MAE <- c(MAE, MAE2)
        return(MAE)
    }

    print("Downloaded published miRNA target interaction data version 7.0")
    print("Currently the servers for miRTarBase are down. hsa or mmu
          data is available here from TimiRGeN.")
return(MAE)
}
