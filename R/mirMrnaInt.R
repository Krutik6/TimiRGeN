#' @title mirMrnaInt
#' @description Create a correlation matrix of all the potential miR-mRNA
#' interactions which could arise between the input microRNA data and and
#' the mRNAs found from the wikiMrna function. The time series DE data will
#' be averaged from the dataframes created by diffExpressRes.
#' @param MAE MultiAssayExperiment which will have the output of
#' mirMrnaInt stored in it. It is recommended to begin a new MAE using
#' MultiAssayExperiment() here so the MAE objects do not get too large.
#' @param miR_express Dataframe from diffExpressRes function of microRNA data.
#' Rownames should be miR gene names and columns should include DE results
#' displaying abundance e.g. log2fc or ave expression. These dataframes should
#' also have gene IDs. This should be stored as an assay within the MAE used
#' in the diffExpressRes function.
#' @param GenesofInterest Input results from wikiMrna. This should be found
#' as an assay within the MAE which was used in the wikiMrna function.
#' @param maxInt Integer. Should be equal to number of samples in both mRNA and
#' miR data e.g. number of different time points. In our example it is 5 because
#' there are 5 time points.
#' @return A large correlation matrix which contains averaged miR-mRNA
#' time series information for every possible miR-mRNA interaction.
#' @export
#' @usage mirMrnaInt(MAE, miR_express, GenesofInterest, maxInt)
#' @examples
#' G <- data.frame(row.names = c("Acaa1a", "Acadm", "Acss1", "Adh1"),
#'                 "D1.Log2FC" = c("-1.2944593","-2.0267432","-2.1934942",
#'                                 "-2.1095853"),
#'                 "D2.Log2FC" = c("-1.1962396","-2.1345451","-1.7699232",
#'                                 "-1.0961674"),
#'                 "D3.Log2FC" = c("0.2738496","-1.9991046","-1.7637549",
#'                                 "-1.6572653"),
#'                 "D7.Log2FC" = c("-0.51765245","-2.20689661","-0.68479699",
#'                                 "-2.06512466"),
#'                 "D14.Log2FC" = c("-0.4510294","-1.1523849","-0.4297012",
#'                                  "-1.1017597"),
#'                 "ID" = c("113868","11364","68738","11522"))
#'
#' MIR <- data.frame(row.names = c("mmu-miR-101a-3p", "mmu-miR-101a-5p",
#'                                 "mmu-miR-101c", "mmu-miR-106a-5p"),
#'                   "D1.Log2FC" = c("-0.0039141722","-0.4328659746",
#'                                   "-0.0038897133", "-0.4161749123"),
#'                   "D2.Log2FC" = c("-0.210605345","-0.600422732",
#'                                   "-0.210574742", "-0.530311376"),
#'                   "D3.Log2FC" = c("-0.315070839","-0.745367163",
#'                                   "-0.315012148", "-0.559274530"),
#'                   "D5.Log2FC" = c("-0.41087763","-0.63952382",
#'                                   "-0.41087876", "-1.03618015"),
#'                   "D14.Log2FC" = c("-0.39466968","-0.60122678",
#'                                    "-0.39461099", "-0.41889698"),
#'                   "ID" = c("387143","387143","100628572","723829"))
#'
#'MAE <- MultiAssayExperiment()
#'
#'MAE <- mirMrnaInt(MAE, miR_express = MIR, GenesofInterest = G,
#'                      maxInt = 5)
mirMrnaInt <- function(MAE, miR_express, GenesofInterest, maxInt){


        if (missing(MAE)) stop('Add MAE object. Output from mirMrnaInt will
                               be added to this MAE as an assay. Please used
                               the diffExpressRes and wikiMrna functions first.')

        if (missing(miR_express)) stop("Add a dataframe with miRNA abundance
                                       values and gene IDs. Please use the
                                       diffExpressRes function on miRNA data.
                                       Output of diffExpressRes will be stored
                                       as an assay within the MAE used in the
                                       diffExpressRes function.")

        if (missing(GenesofInterest)) stop("Add a dataframe containing mRNAs
                                           found in the wikipathway of interest
                                           and in the input mRNA data. Please
                                           use the wikiMrna function. Results
                                           from wikiMrna will be stored as
                                           an assay the MAE used in the
                                           wikiMrna function.")

        if (missing(maxInt)) stop("Integer. Number of time points your
                                  time course goes along.")

        #Make matrix of the names of miRs and mRNAs of interest
        allnames <- expand.grid(rownames(miR_express),
                                rownames(GenesofInterest),
                                stringsAsFactors = FALSE)

        # Apply CorrTable to this
        interactions_df <- do.call(rbind, apply(allnames, 1,
                                                corrTable,
                                                miR_express,
                                                GenesofInterest,
                                                maxInt = maxInt))

        MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                      interactions_df = interactions_df)))

        MAE <- c(MAE, MAE2)
return(MAE)
}
