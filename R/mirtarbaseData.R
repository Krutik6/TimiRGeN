#' @title mirtarbaseData
#' @description Produces a list of mRNA-miR interactions with functional
#' validation from miRTarBase (weaker validation methods are not included).
#' @param mirtarbase Output from dloadMirtarbase. Download from here
#' http://mirtarbase.mbc.nctu.edu.tw/php/download.php.
#' @param species Species being studies e.g. hsa or mmu.
#' @param MAE Multiassayexperiment object.
#' @return Dataframe of mRNA-miR interactions with functional evidence.
#' @export
#' @importFrom stringr str_detect
#' @importFrom dplyr filter
#' @import tidyverse
#' @usage mirtarbaseData(MAE, mirtarbase, species = '')
#' @examples
#' MAE <- MultiAssayExperiment(list("mirtarbase" = data.frame(
#'                             miRTarBase.ID = c("MIRT006298", "MIRT737375"),
#'                             miRNA = c("hsa-miR-200b-3p", "hsa-miR-1321"),
#'                             Species..miRNA. = c("Homo sapiens",
#'                                                 "Homo sapiens"),
#'                             Target.Gene = c("E2F3", "STX1A"),
#'                             Target.Gene..Entrez.ID.= c(1871, 6804),
#'                             Species..Target.Gene.=c("Homo sapiens",
#'                                                      "Homo sapiens"),
#'                            Experiments = c(" Flow//Luciferase reporter assay//Western blot",
#'                                            "PAR-CLIP"),
#'                            Support.Type = c("Functional MTI",
#'                                             "Functional MTI (Weak)"),
#'                            References..PMID. = c(22144583, 26701625))))
#'
#' MAE <- mirtarbaseData(MAE = MAE, mirtarbase = assay(MAE, 1),
#'                       species = "hsa")
mirtarbaseData <- function(MAE, mirtarbase, species){

  if (missing(MAE)) stop('Add MAE object.')

  if (missing(mirtarbase)) stop('Input mirtarbase file.')

  if (missing(species)) stop('e.g hsa or mmu.')

  miRNA <- NULL

  # Extract species specific data
  miRTarBase_s <- mirtarbase %>% dplyr::filter(stringr::str_detect(
                                                                miRNA, species))

  # Remove weak functional MTI
  miRTarBase_Fun <- miRTarBase_s[which(
                               miRTarBase_s$Support.Type == 'Functional MTI'),]


  # Organise output in a structured way
  miRTarBase_df <- data.frame(miRTarBase_Interactions = paste(
                                miRTarBase_Fun$miRNA, ':',
                                miRTarBase_Fun$Target.Gene, sep = ''),
                                miRTarBase_microRNA = miRTarBase_Fun$miRNA,
                                miRTarBase_mRNA = miRTarBase_Fun$Target.Gene)

  MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                             'miRTarBase_res' = miRTarBase_df)))

  MAE <- c(MAE, MAE2)

return(MAE)
}
