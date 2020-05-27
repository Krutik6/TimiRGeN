#' @title genesList
#' @description Will produce a list of nested dataframes. For combined
#' analysis method = c, and for separated analysis method = s. Combined
#' analysis colnames should be 'timepoint.resulttype'.
#' genesList will make new dataframes separated by 'timepoint.'.
#' Separated analysis colnames should be 'genetype_timepoint.resulttype'.
#' genesList will make separate lists for 'genetype_'. Each of these lists
#' will have dataframes which have been made by separating by 'timepoint.'.
#' Make sure to follow colname nomenclature carefully.
#' @param MAE Multi Assay Experiment object.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param genetic_data If 'c', this should be a dataframe with miR and mRNA
#' information together.
#' @param timeString If 'c', this should be a commmon string representing
#' 'timepoint' e.g. for H.1, H.10, H.20, timeString = 'H'.
#' @param miR_data If 's', a dataframe of microRNA data. Rownames are genes
#' and colnames are: genetype_timepoint.resulttype or timepoint.resulttype.
#' timepoint.resulttype should be the same in mRNA and miR data.
#' @param mRNA_data If 's', a dataframe of mRNA data. Rownames are genes
#' and colnames are: genetype_timepoint.resulttype or timepoint.resulttype.
#' @return A list of dataframes separated by features in the column names which
#' can be stored in the metadata area of an MAE object.
#' @export
#' @importFrom stringr str_extract
#' @usage genesList(MAE, method, genetic_data, timeString, miR_data, mRNA_data)
#' @examples
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
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
#' MAE <- combineGenes(MAE, miR_data = assay(MAE, 1),
#'                     mRNA_data = assay(MAE, 2))
#'
#' MAE <- genesList(MAE, method = 'c', genetic_data = assay(MAE, 5),
#'                  timeString = 'D')
genesList <- function(MAE, method, genetic_data, timeString, miR_data,
                      mRNA_data){

    if (missing(MAE)) stop('Add MAE a object')

    if (method == 'c') {

        if (missing(genetic_data)) stop('Input combined miR and mRNA data.
                                        Colnames structure should be
                                        timepoint.resulttype.')

        if (missing(timeString)) stop('Input timepoint. E.g. if colnames were
                                      D1, D2, D3; then timeString = D.')

        metadata <- `metadata<-` <- NULL

        G <- genetic_data

        # Change random string into common string
        colnames(G) <- gsub(x = colnames(G), pattern = timeString,
                            replacement = 'TP')

        # Split into lists by TP
        X <- lapply(split.default(G, sub("(TP\\d+).*", "\\1",
                                         names(G))), as.list)

        # Relplace TP with timeString
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
        if (missing(miR_data)) stop('Input miR data. Colnames structure
                                    should be genetype_timepoint.resulttype.')

        if (missing(mRNA_data)) stop('Input mRNA data. Colnames structure
                                     should be genetype_timepoint.resulttype.')


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

    } else{ stop('Please insert method c for combined analysis or s for
                 seperate analysis')}
}
