#' @title genesList
#' @description Produces a list of nested dataframes. The list will depend on
#' the type of analysis that is to be conducted.
#' For combined analysis method = "c", and for separated analysis method = "s".
#'
#' In combined analysis colnames should be 'timepoint.resulttype'.
#' genesList will make new dataframes separated at 'timepoint.'.
#'
#' In separated analysis colnames should be 'genetype_timepoint.resulttype'.
#' genesList will make separate lists for each 'genetype_', and these lists
#' will have dataframes which have been made by separating at 'timepoint.'.
#'
#' Make sure to follow colname nomenclature carefully. Please refer to
#' the vignette for more details on the nomenclature.
#' @param MAE MultiAssayExperiment which will store the output from genesList.
#' It is recommended to use the MAE which stores output from combineGenes
#' (combined analysis) or addPrefix (separated analysis).
#' @param method Either "c" or "s", respectively for combined or separated
#' analysis.
#' @param genetic_data If "c", this should be a dataframe with miR and mRNA
#' information together. This is the output from the combineGenes function
#' and will be stored as an assay within the MAE used in the combineGenes
#' function.
#' @param timeString If "c", this should be a common string representing
#' 'timepoints' e.g. for H.1, H.10, H.20, timeString = "H".
#' @param miR_data If "s", a dataframe of microRNA data. Rownames are genes
#' and colnames are: genetype_timepoint.resulttype. Column names should be
#' the same in mRNA and miR data. miR_data is from the addPrefix function,
#' and will be stored as an assay within the MAE used in addPrefix.
#' @param mRNA_data If "s", a dataframe of mRNA data. Rownames are genes
#' and colnames are: genetype_timepoint.resulttype. Column names should be the
#' same in mRNA and miR data. mRNA_data is from the addPrefix function,
#' and will be stored as an assay within the MAE used in addPrefix.
#' @return A list of dataframes separated by features in the column names.
#' Output will be stored as metadata in the input MAE.
#' @export
#' @importFrom stringr str_extract
#' @usage genesList(MAE, method, genetic_data, timeString, miR_data, mRNA_data)
#' @examples
#' miR <- mm_miR[1:50,]
#'
#' mRNA <- mm_mRNA[1:100,]
#'
#' MAE <- startObject(miR = mm_miR, mRNA = mm_mRNA)
#'
#'# For separated analysis
#'
#' MAE <- addPrefix(MAE = MAE, gene_df = assay(MAE, 1),
#'                  prefixString = "miR")
#'
#' MAE <- addPrefix(MAE = MAE, gene_df = assay(MAE, 2),
#'                  prefixString = "mRNA")
#'
#' MAE <- genesList(MAE, method = "s", miR_data = assay(MAE, 3),
#'                  mRNA_data = assay(MAE, 4))
#'
#'# For combined analysis
#'
#' MAE <- combineGenes(MAE, miR_data = assay(MAE, 1),
#'                     mRNA_data = assay(MAE, 2))
#'
#' MAE <- genesList(MAE, method = 'c', genetic_data = assay(MAE, 3),
#'                  timeString = 'D')
genesList <- function(MAE, method, genetic_data, timeString, miR_data,
                      mRNA_data){

    if (missing(MAE)) stop('MAE is missing. Add MAE. Please use combineGenes or addPrefix first.')

    if (method == 'c') {


        if (missing(genetic_data)) stop('genetic_data is missing. Input combined miR and mRNA data. Colnames structure should be timepoint.resulttype. Please use the combineGenes function first. Output of combineGenes should be stored as an assay within the MAE used in the combineGenes function.')

        if (missing(timeString)) stop('timeString is missing. Input unit of time. E.g. if colnames were D1.log2fc, D2.log2fc, D3.log2fc; then timeString = "D".')

        metadata <- `metadata<-` <- NULL

        G <- genetic_data

        # Change random string into common string
        colnames(G) <- gsub(x = colnames(G), pattern = timeString,
                            replacement = 'TP')

        # Split into lists by TP
        X <- lapply(split.default(G, sub("(TP\\d+).*", "\\1",
                                         names(G))), as.list)

        # Replace TP with timeString
        names(X) <- gsub(names(X), pattern = 'TP', replacement = timeString)

        # Separate into nested data frames revolving around only the time point
        L1 <- lapply(X, data.frame, stringsAsFactors = FALSE)

        L2 <- lapply(L1, function(DF) {rownames(DF) <- rownames(
                                                            genetic_data); DF})
        L3 <- L2[gtools::mixedsort(names(L2))]

        geneslist <- lapply(L3, function(x) {
            colnames(x) <- sub(x = colnames(x),
                               pattern = 'TP',
                               replacement = timeString)
            x}
        )

        # Store into MAE object
        metadata(MAE)[["geneslist"]] <- geneslist

        return(MAE)

    } else if (method == 's') {
        if (missing(miR_data)) stop('miR_data is missing. Input miR data with prefixes. Colnames structure should be genetype_timepoint.resulttype. Please use the addPrefix function first. Output of the addPrefix function should be stored as an assay within the MAE used in the addPrefix function.')

        if (missing(mRNA_data)) stop('mRNA_data is missing. Input mRNA data with prefixes. Colnames structure should be genetype_timepoint.resulttype. Please use the addPrefix function first. Output of the addPrefix function should be stored as an assay within the MAE used in the addPrefix function.')


        # get data from MAE object
        miR_data <- as.data.frame(miR_data)

        mRNA_data <- as.data.frame(mRNA_data)

        # Combine this data into a single list
        genedata <- list(miR_data = miR_data, mRNA_data = mRNA_data)

        # Split into lists by data type and time point
        geneslist <- lapply(genedata, function(df)
        {
            Unigenes <- unique(stringr::str_extract(names(df),"\\S+\\."))
            List <- lapply(Unigenes,function(name){return(df[,grep(
                                                                   name,
                                                                   names(df),
                                                                   fixed=TRUE
                                                               )])})
            names(List) <- Unigenes

            return(List)
        })

        # Store into MAE object
        metadata(MAE)[["geneslist"]] <- geneslist

        return(MAE)

    } else print('Enter c or s as method.')
}
