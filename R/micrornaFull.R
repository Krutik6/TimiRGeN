#' @title micrornaFull
#' @description Changes many potential miR naming systems into miR names which
#' adhere to the TimiRGeN friendly miR naming system in one step.
#' This is an internal function for getIDsMir functions.
#' @param miRdf miR dataframe with genes as rownames and results from DE as
#' column names. Will be the miR input for the getIdsMir function used.
#' @param species Species e.g. "hsa" or "mmu".
#' @return hsa-miR-140-5p would be changed to microRNA 140-1
#' @noRd
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
