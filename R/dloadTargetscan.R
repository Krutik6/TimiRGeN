#' @title dloadTargetscan
#' @description Downloads most recent version (7.2) of predicted targets from
#' the targetscan database
#' http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi.
#' miR-mRNA interactions from the species of interest will be extracted.
#' @param MAE MultiAssayExperiment which will store the downloaded targetscan
#' data. It is recommended to use the MAE which was used in
#' the mirMrnaInt function.
#' @param species Species of interest e.g "hsa" or "mmu."
#' @return Dataframe of species specific predicted mRNA-miR interactions. Output
#' will be stored as an assay in the input MAE.
#' @export
#' @usage dloadTargetscan(MAE, species)
#' @examples
#' \dontrun{
#' MAE <- MultiAssayExperiment()
#'
#' MAE <-dloadTargetscan(MAE, "mmu")
#' }
dloadTargetscan <- function(MAE, species){

  if (missing(MAE)) stop('MAE is missing. Add MultiAssayExperiment which will have downloaded data stored in it. Please use the mirMrnaInt function first.')

  if (missing(species)) stop('species is missing. Add initials of the species of interest e.g "hsa" or "mmu".')

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

  Transcript <- NULL

  TargetScans <- as.data.frame(TargetScans)

  # Only keep selected species microRNAs
  TargetScans_s <- TargetScans %>% dplyr::filter(stringr::str_detect(
                                                  Transcript, species))


  if (species == 'hsa') {
    # Specify how to organise data
    Targetscans_df <- data.frame(Targetscans_Interactions =
                                   paste(TargetScans_s$Transcript,
                                         ':',
                                         TargetScans_s$ID,sep = ''),
                                 Targetscans_miR = TargetScans_s$Transcript,
                                 Targetscans_mRNA = TargetScans_s$ID)

    MAE2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(
      'Targetscans_res' = Targetscans_df))))

    MAE <- suppressWarnings(c(MAE, MAE2))

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

    #Store in a MAE
    MAE2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(
                                         'Targetscans_res' = Targetscans_df))))

    MAE <- suppressWarnings(c(MAE, MAE2))

    return(MAE)

  }
}

