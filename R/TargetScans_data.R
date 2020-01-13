#' @title TargetScans_data
#' @description Produces a list of mRNA-miR predicted interactions from miRDB
#'for the species of interest from targetscans. Requires the tidyverse
#'package.
#' @param targetScan File of non-conserved miR-mRNA predicted targets.
#' Downloaded from here
#'http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi.
#' @param species Species we are interested in e.g hsa or mmu.
#' @return List of predicted mRNA-miR interactions and associated information.
#' @export
#' @usage TargetScans_data(targetScan, species = '')
#' @examples
#'library(tidyverse)
#'data.frame(row.names = c("215", "216", "286", "287"),
#'"Gene" = c("ENSG00000148584.10", "ENSG00000148584.10",
#'"ENSG00000109576.9", "ENSG00000109576.9"),
#'"ID" =c("A1CF", "A1CF", "AADAT", "AADAT"),
#'"Gene.1" = c("ENST00000374001.2", "ENST00000374001.2",
#'"ENST00000337664.4", "ENST00000337664.4"),
#'"Symbol" = c("10090", "10090", "10090", "10090"),
#'"Transcript" = c("mmu-miR-429-3p", "mmu-miR-374c-5p",
#'"mmu-miR-182-5p",
#'"mmu-miR-340-5p")) -> TargetScans
#'TargetScans_data(targetScan = TargetScans,
#'species = 'mmu') -> TargetScans_mmu
TargetScans_data <- function(targetScan, species){
if (missing(targetScan)) stop('Input targetscans file.');
if (missing(species)) stop('e.g hsa or mmu.');
Transcript <- NULL
targetScan %>%
filter(str_detect(Transcript, species)) -> TargetScans_s
if (species == 'hsa') {
Targetscans_df <- data.frame(Targetscans_Interactions =
paste(TargetScans_s$Transcript,
':', TargetScans_s$ID,sep = ''),
Targetscans_miR = TargetScans_s$Transcript,
Targetscans_mRNA = TargetScans_s$ID)
return(Targetscans_df)
}else if (species == 'mmu') {
tolower(TargetScans_s$ID) -> TargetScans_s$ID2
firstup <- function(x) {
substr(x, 1, 1) <- toupper(substr(x, 1, 1))
x
}
firstup(TargetScans_s$ID2) -> TargetScans_s$ID2
Targetscans_df <- data.frame(Targetscans_Interactions =
paste(TargetScans_s$Transcript,
':', TargetScans_s$ID2,sep = ''),
Targetscans_miR = TargetScans_s$Transcript,
Targetscans_mRNA = TargetScans_s$ID2)
return(Targetscans_df)
} else (stop('Add hsa or mmu as species.'))
}
