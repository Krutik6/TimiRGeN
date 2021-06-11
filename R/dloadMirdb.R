#' @title dloadMirdb
#' @description Downloads most recent version (6.0) of predicted targets from
#' the mirdb database http://mirdb.org/download.html. This will take some time.
#' miR-mRNA interactions from the species of interest will be extracted.
#' Species of interests associated org.db package must be loaded beforehand.
#' @param MAE MultiAssayExperiment which will store downloaded mirDB data.
#' It is recommended to use the MAE which was used in the mirMrnaInt function.
#' @param species Species of interest e.g. "hsa" or "mmu".
#' @param orgDB Organism db package specific for species of interest.
#' @return A dataframe of predicted, species specific mRNA-miR interactions.
#' Will be stored as an assay in the input MAE.
#' @export
#' @usage dloadMirdb(MAE, species, orgDB)
#' @examples
#' \dontrun{
#'
#' library(org.Mm.eg.db)
#'
#' MAE <- MultiAssayExperiment()
#'
#' MAE <-dloadMirdb(MAE, 'mmu', org.Mm.eg.db)
#'
#' }
dloadMirdb <- function(MAE, species, orgDB){

  if (missing(MAE)) stop('MAE is missing. Add MAE which will have output from dloadMirdb stored within. Please use mirMrnaInt first.')

  if (missing(species)) stop('species is missing. Add initials of the species of interest e.g "hsa" or "mmu".')

  if (missing(orgDB)) stop('orgDB is missing. Add a species library e.g. org.Hs.eg.db.')

  miR <- NULL

  # download targetscans data
  download.file("http://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
                'MIRDB.gz')

  # read file from directory
  miRDB <- read.table(gzfile("MIRDB.gz"), header = FALSE, fill = TRUE)

  # remove file
  file.remove("MIRDB.gz")

  names(miRDB) <- c('miR', 'Target', 'Score')

  # Isolate mouse/ human specific miRs
  miRDB_s <- miRDB %>% dplyr::filter(stringr::str_detect(miR, species))

  # Get human readable gene names
  miRDB_mRNA <- suppressMessages(suppressWarnings(
                                  clusterProfiler::bitr(miRDB_s$Target,
                                                       fromType = 'REFSEQ',
                                                       toType = 'SYMBOL',
                                                       OrgDb = orgDB)))

  # Link species specific miRs to their targets
  mirDB_merged <- merge(miRDB_s, miRDB_mRNA,
                        by.x = 'Target',
                        by.y = 'REFSEQ',
                        all = TRUE)

  # Organise this information
  miRDB_df <- data.frame(miRDB_Interactions = paste(mirDB_merged$miR,
                                                    ':',
                                                    mirDB_merged$SYMBOL,
                                                    sep = ''),
                         miRDB_miR = mirDB_merged$miR,
                         miRDB_mRNA = mirDB_merged$SYMBOL)


  #Store in MAE
  MAE2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list('miRDB_res' = miRDB_df)
  )))

  MAE <- suppressWarnings(c(MAE, MAE2))

  return(MAE)
}

