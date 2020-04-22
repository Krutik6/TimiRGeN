#' @title micrornaFull
#' @description Changes miR names into microRNA names in one step.
#' This has less flexibility that MicroRNA_partial.
#' @param miRdf miR dataframe.
#' @param species Species e.g. hsa or mmu.
#' @return hsa-miR-140-5p would be changed to microRNA 140-1
#' @export
#' @examples
#' miR <- mm_miR
#' microRNA_genenames <- micrornaFull(rownames(miR), species = "mmu")
micrornaFull <- function(miRdf, species){

    if (missing(miRdf)) stop('Input rowname or column with microRNA
                              names.');
    if (missing(species)) stop('e.g hsa or mmu.');

    x <- gsub(miRdf, pattern = paste0(species,'-miR-'),
              replacement = 'microRNA ')

    x <- gsub(x, pattern = species, replacement = 'microRNA ')

    x <- gsub(x, pattern = '5p', replacement = '1')

    x <- gsub(x, pattern = '3p', replacement = '2')

    x <- gsub(x, pattern = '-let', replacement = 'let')

return(x)
}
