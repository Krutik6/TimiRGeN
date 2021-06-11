#' @title dloadMirtarbase
#' @description Downloads most recent version (8.0) of functional targets from
#' the miRTarBase database http://mirtarbase.cuhk.edu.cn/php/download.php
#' . Species specific miR-mRNA interactions which do not have 'weak' evidence
#' are used.
#' @param MAE MultiAssayExperiment which will store the downloaded mirtarbase
#' data. It is recommended to use the MAE which was used in
#' the mirMrnaInt function.
#' @param species Species of interest e.g. "hsa" or "mmu".
#' @return Dataframe of species specific miR-mRNA interactions with
#' strong functional evidence. Output will be stored as an assay in the input
#' MAE.
#' @export
#' @usage dloadMirtarbase(MAE, species)
#' @importFrom readxl read_excel
#' @examples
#' MAE <- MultiAssayExperiment()
#'
#' MAE <- dloadMirtarbase(MAE, "mmu")
#'
dloadMirtarbase <- function(MAE, species){

  if (missing(MAE)) stop('MAE is missing. Add MultiAssayExperiment so data from dloadMirtarbase can be stored. Please use the mirMrnaInt function first.')

  if (missing(species)) stop('species is missing. Add initials of the species of interest e.g "hsa" or "mmu."')


  miRTarBase <- miRNA <- NULL

  # Retreive data
  miRTarBase <- TimiRGeN::miRTarBase

  # Extract species specific data
  miRTarBase_s <- miRTarBase %>% dplyr::filter(stringr::str_detect(miRNA, species))

  # Organise output in a structured way
  miRTarBase_df <- data.frame(miRTarBase_Interactions = paste(
                              miRTarBase_s$miRNA,
                              ':',
                              miRTarBase_s$Target.Gene, sep = ''),
                              miRTarBase_microRNA = miRTarBase_s$miRNA,
                              miRTarBase_mRNA = miRTarBase_s$Target.Gene)

  #Store in a MAE
  MAE2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(
                              'miRTarBase_res' = miRTarBase_df))))

  MAE <- suppressWarnings(c(MAE, MAE2))

  return(MAE)
}
