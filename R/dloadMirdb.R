#' @title dloadMirdb
#' @description Downloads most recent version (6.0) of predicted targets from
#' mirdb database http://mirdb.org/download.html. This will take some time.
#' miR-mRNA interactions from the species of interest will be extracted and
#' stored in a MAE object as an assay. Species associated org.db package must
#' be loaded.
#' @param MAE MultiAssayExperiment which will have downloaded mirDB data stored
#' within it. It is recommended to use the MAE which was created using the
#' mirMrnaInt function.
#' @param species Species of interest e.g. "hsa" or "mmu".
#' @param orgDB organism db specific for species being looked into. Specific
#' db must be loaded before use.
#' @return A dataframe of species specific mRNA-miR interactions, and other
#' information associated with these interactions.
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

  if (missing(MAE)) stop('Add MAE object which will have output from
                          dloadMirdb stored in. Please use mirMrnaInt first.')

  if (missing(species)) stop('initials of the species of interest
                               e.g hsa or mmu.')

  if (missing(orgDB)) stop('Add a speices library e.g. org.Hs.eg.db.')

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
  miRDB_mRNA <- suppressMessages(clusterProfiler::bitr(miRDB_s$Target,
                                                       fromType = 'REFSEQ',
                                                       toType = 'SYMBOL',
                                                       OrgDb = orgDB))

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
  MAE2 <- suppressMessages(MultiAssayExperiment(list('miRDB_res' = miRDB_df)
  ))

  MAE <- c(MAE, MAE2)

  return(MAE)
}

