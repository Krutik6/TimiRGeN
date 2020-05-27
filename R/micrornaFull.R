#' @title micrornaFull
#' @description Changes miR names into microRNA names in one step. Internal
#' function for getIDsMir functions.
#' @param miRdf miR dataframe.
#' @param species Species e.g. hsa or mmu.
#' @return hsa-miR-140-5p would be changed to microRNA 140-1
micrornaFull <- function(miRdf, species){

    #  microRNA > miR
    x <- gsub(miRdf, pattern = paste0(species,'-miR-'),
              replacement = 'microRNA ')

    # microRNA > hsa/mmu
    x <- gsub(x, pattern = species, replacement = 'microRNA ')

    # replace 5p or 3p to 1 and 2 respectfully
    x <- gsub(x, pattern = '5p', replacement = '1')

    x <- gsub(x, pattern = '3p', replacement = '2')

    # add - infront of let genes
    x <- gsub(x, pattern = '-let', replacement = 'let')

return(x)
}
