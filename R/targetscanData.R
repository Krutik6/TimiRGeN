#' @title targetscanData
#' @description Produces a list of mRNA-miR predicted interactions from
#' dloadTargetscans, for the species of interest.
#' @param MAE MultiAssayExperiments.
#' @param targetScan File of non-conserved miR-mRNA predicted targets.
#' Downloaded from here
#'http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi,
#'using dloadTargetscans.
#' @param species Species we are interested in e.g hsa or mmu.
#' @return List of predicted mRNA-miR interactions and associated information.
#' @export
#' @importFrom stringr str_detect
#' @importFrom clusterProfiler bitr
#' @import tidyverse
#' @usage targetscanData(MAE, targetScan, species = '')
#' @examples
#'MAE <- MultiAssayExperiment(list(TargetScans = data.frame(
#'                          row.names = c("215", "216", "286", "287"),
#'                         "Gene" = c("ENSG00000148584.10",
#'                                    "ENSG00000148584.10",
#'                                    "ENSG00000109576.9",
#'                                    "ENSG00000109576.9"),
#'                          "ID" =c("A1CF", "A1CF", "AADAT", "AADAT"),
#'                          "Gene.1" = c("ENST00000374001.2",
#'                                       "ENST00000374001.2",
#'                                       "ENST00000337664.4",
#'                                       "ENST00000337664.4"),
#'                         "Symbol" = c("10090", "10090", "10090", "10090"),
#'                        "Transcript" = c("mmu-miR-429-3p", "mmu-miR-374c-5p",
#'                                     "mmu-miR-182-5p","mmu-miR-340-5p"))))
#'
#'TargetScans_df <- targetscanData(MAE, targetScan = assay(MAE, 1),
#'                                   species = 'mmu')
targetscanData <- function(MAE, targetScan, species){

    if (missing(MAE)) stop('Use MultiAssayExperiment.')

    if (missing(targetScan)) stop('Output from dloadTargetscans.')

    if (missing(species)) stop('e.g hsa or mmu.')

    Transcript <- NULL

    targetScan <- as.data.frame(targetScan)

    # Only keep selected species microRNAs
    TargetScans_s <- targetScan %>% dplyr::filter(stringr::str_detect(
                                                           Transcript, species))

    if (species == 'hsa') {
        # Specify how to organise data
        Targetscans_df <- data.frame(Targetscans_Interactions =
                                     paste(TargetScans_s$Transcript,
                                     ':',
                                     TargetScans_s$ID,sep = ''),
                                     Targetscans_miR = TargetScans_s$Transcript,
                                     Targetscans_mRNA = TargetScans_s$ID)

        MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                           'Targetscans_res' = Targetscans_df)))

        MAE <- c(MAE, MAE2)

    return(MAE)

    }else if (species == 'mmu') {

        # Change mRNA names to capitals
        TargetScans_s$ID2 <- tolower(TargetScans_s$ID)

        firstup <- function(x) {
                   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
                    x
        }
        TargetScans_s$ID2 <- firstup(TargetScans_s$ID2)

        # Specify how to organise data
        Targetscans_df <- data.frame(Targetscans_Interactions =
                                    paste(TargetScans_s$Transcript,
                                    ':',
                                    TargetScans_s$ID2,sep = ''),
                                    Targetscans_miR = TargetScans_s$Transcript,
                                    Targetscans_mRNA = TargetScans_s$ID2)

        MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                        'Targetscans_res' = Targetscans_df)))
        MAE <- c(MAE, MAE2)
    return(MAE)

} else (stop('Add hsa or mmu as species.'))
}
