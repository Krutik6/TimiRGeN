#' @title CreateClusters
#' @description Create soft clusters to assess how your datasets change in
#' gene abundance during your time course in different wikipathways.
#' This function is to be used before ClusterCheck and Quickfuzz.
#' @param method Either C or S for combined or separate analysis.
#' @param Percent_matrix Output after TurnPercent.
#' @param no.clusters Number of clusters to create, the default is 5.
#' @param Data_string Only for use in S analysis. Insert the prefix string e.g.
#' mRNA or miR.
#' @param Variance Numeric from 0-1 to control strictness of filtering. Higher
#' Variance means wikipathways that are undergoing a higher change in the
#' data will be clustered.
#' @importFrom Mfuzz filter.std standardise mestimate mfuzz
#' @return Clusters: A list to be used as the input in plot functions.
#' Mfuzzdata: An input for QuickFuzz. Mfuzzdata: An ExpressionSet object
#' to be input for Mfuzz clustering.
#' @export
#' @usage CreateClusters(method, Percent_matrix, no.clusters,
#' Data_string = '', Variance)
#' @examples
#' library(Mfuzz)
#' w <- list(WP1 = c("1071", "11303", "11806", "11808", "11812", "11813",
#'"11814", "11816", "13122", "13350", "15357", "15450",
#'"16816",  "16835", "16956", "16971", "17777", "18830",
#'"20652", "20778"),
#'WP10 = c("11651", "12702", "14784", "16186", "16198", "16199",
#'"16367", "16451", "16453", "18708", "19247", "20416", "20846",
#'"20848", "20850", "20851", "26395", "26396", "26413", "26417",
#'"269523", "384783", "54721", "81601"),
#'WP103 = c("110196", "13121",  "13360", "14137", "15357", "16987",
#'"17855", "18194", "192156", "20775", "208715", "235293",
#'"319554", "66234", "68603"),
#'WP108 = c("107585",  "107869", "109079", "109815", "114679", "12916",
#'"13370", "13371",  "14080", "14281", "14775", "14776",
#'"14778", "16476", "18024",  "18033", "18986", "19697",
#'"19946", "20226", "20341", "20342",  "20363", "20364",
#'"20683", "20687", "20768", "211006", "214580",
#'"223776", "232223", "26462", "27361", "28042", "280621",
#'"433003", "50493", "50880", "619883", "625249", "625281",
#'"65967", "664969",  "69227", "71787", "71984", "72657",
#'"74777", "75420", "75512",  "80795", "9403"),
#'WP113 = c("12159", "12387", "12393", "12399"))
#'e <- list(D1 = c("1071", "11303", "11806", "11808", "223776", "232223",
#'"26462", "27361", "17855", "18194", "192156" ),
#'D2 = c("19946", "20226", "20341", "20342",  "20363", "20364",
#'"20683", "20687", "20768", "211006", "214580", "16971",
#'"17777", "18830"),
#'D3 = c("1071", "11303", "11806", "11808", "11812", "11813",
#'"11814", "11816", "13122", "13350", "20652", "20778"))
#'WikiMatrix(e_list = e, wp_list = w) -> wmat
#'TurnPercent(wikiMatrix = wmat, rowInt = 4) -> pmat
#'CreateClusters(method = 'c', Percent_matrix = pmat,
#'no.clusters = 3, Variance = 0)
CreateClusters <- function(method, Percent_matrix, no.clusters = 5, Data_string,
        Variance = 0){
        as.data.frame(t(Percent_matrix)) -> df
        df$Total <- NULL
        df[vapply(df, is.factor, logical(1))] <- lapply(df[vapply(
        df, is.factor, logical(1))], function(x) as.numeric(as.character(x)))
        round(df, 0) -> df
        na.omit(df) -> df
        if (method == 's') {
        df[, grepl(Data_string, names(df))] -> df2
        } else if (method == 'c') {
        df -> df2
        } else print('Select s for seperated analysis or c for combined
        analysis')
        Eset <- new('ExpressionSet', exprs = as.matrix(df2))
        Eset_sd <- filter.std(Eset, min.std = Variance)
        Eset_st <- standardise(Eset_sd)
        m <- mestimate(Eset_st)
        cl <- mfuzz(Eset_st, centers = no.clusters, m=m)
        cl$membership -> ClusterData
        list(ClusterData = ClusterData,
        Clusters = cl,
        Mfuzzdata = Eset_st) -> clusterlist
return(list2env(clusterlist, .GlobalEnv))
}
