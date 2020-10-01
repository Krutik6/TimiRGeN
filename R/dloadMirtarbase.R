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
#' \dontrun{
#'
#' MAE <- MultiAssayExperiment()
#'
#' MAE <- dloadMirtarbase(MAE, "mmu")
#'
#' }
dloadMirtarbase <- function(MAE, species){

  if (missing(MAE)) stop('
                         MAE is missing.
                         Add MultiAssayExperiment so data from dloadMirtarbase
                         can be stored. Please use the mirMrnaInt function
                         first.')

  if (missing(species)) stop('
                             species is missing.
                             Add initials of the species of interest
                             e.g "hsa" or "mmu."')

  # download miRTarBase data
  download.file(
    'http://mirtarbase.cuhk.edu.cn/cache/download/8.0/miRTarBase_MTI.xlsx',
    'miRTarBase.xlsx')

  # read using readxl
  miRTarBase <- readxl::read_excel('miRTarBase.xlsx')

  # write it as a CSV file so it can be imported back in the intended structure
  write.csv(miRTarBase, 'miRTarBase.csv')

  miRTarBase <- read.csv('miRTarBase.csv', row.names = 1)

  # Remove multiple files
  file.remove("miRTarBase.xlsx")

  file.remove("miRTarBase.csv")

  miRNA <- NULL

  # Extract species specific data
  miRTarBase_s <- miRTarBase %>% dplyr::filter(stringr::str_detect(
                                               miRNA, species))

  # Remove weak functional MTI
  miRTarBase_Fun <- miRTarBase_s[which(
                               miRTarBase_s$Support.Type == 'Functional MTI'),]


  # Organise output in a structured way
  miRTarBase_df <- data.frame(miRTarBase_Interactions = paste(
                                        miRTarBase_Fun$miRNA,
                                        ':',
                                        miRTarBase_Fun$Target.Gene, sep = ''),
                              miRTarBase_microRNA = miRTarBase_Fun$miRNA,
                              miRTarBase_mRNA = miRTarBase_Fun$Target.Gene)

  #Store in a MAE
  MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                            'miRTarBase_res' = miRTarBase_df)))

  MAE <- c(MAE, MAE2)

  return(MAE)
}


