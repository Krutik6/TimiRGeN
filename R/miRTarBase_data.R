#' @title miRTarBase_data
#' @description Produces a list of mRNA-miR interactions with functional
#' validation (weaker validation methods are not included) from miRDB for the
#'species of interest from targetscans.
#' @param mirtarbase Datafile of the most recent experimental mRNA-miR
#' interactions. Download from here
#' http://mirtarbase.mbc.nctu.edu.tw/php/download.php.
#' @param species Species being studies e.g. hsa or mmu.
#' @return Dataframe of mRNA-miR interactions with functional evidence.
#' @export
#' @usage miRTarBase_data(mirtarbase, species = '')
#' @examples
#' library(tidyverse)
#' data.frame(row.names = c("506668", "506669", "506670"),
#' "miRTarBase.ID" = c("MIRT000005", "MIRT000005", "MIRT000016"),
#' "miRNA" = c("mmu-miR-124-3p", "mmu-miR-124-3p", "mmu-miR-210-3p"),
#' "Species..miRNA." = c(rep("Mus musculus", 3)),
#' "Target.Gene" = c("Itgb1", "Itgb1", "Tcf7l2"),
#' "Target.Gene..Entrez.ID." = c("16412", "16412", "21416"),
#' "Species..Target.Gene." = c(rep("Mus musculus", 3)),
#' "Experiments" = c("Luciferase reporter assay//Microarray//qRT-PCR",
#'"Luciferase reporter assay//qRT-PCR//Western blot//Reporter assay;Microarray",
#'"Luciferase reporter assay//Western blot//Microarray"),
#'"Support.Type" = c(rep("Functional MTI", 3)),
#'"References..PMID." = c("18042700", "18619591",
#'"20492721")) -> miRTarBase
#'miRTarBase_data(mirtarbase = miRTarBase, species = 'mmu') -> miRTarBase_mmu
miRTarBase_data <- function(mirtarbase, species){
        if (missing(mirtarbase)) stop('Input targetscans file.');
        if (missing(species)) stop('e.g hsa or mmu.');
        miRNA <- NULL
        mirtarbase %>%
        filter(str_detect(miRNA, species)) -> miRTarBase_s
        miRTarBase_s[which(miRTarBase_s$Support.Type == 'Functional MTI'),
        ] -> miRTarBase_Fun
        miRTarBase_df <- data.frame(miRTarBase_Interactions =
        paste(miRTarBase_Fun$miRNA, ':',
        miRTarBase_Fun$Target.Gene, sep = ''),
        miRTarBase_microRNA = miRTarBase_Fun$miRNA,
        miRTarBase_mRNA = miRTarBase_Fun$Target.Gene)
return(miRTarBase_df)
}
