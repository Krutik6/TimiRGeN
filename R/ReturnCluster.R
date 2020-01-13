#' @title ReturnCluster
#' @description Will retrive information about which wikipathways fitted best
#' with a specific cluster.
#' @param ClusterData Output after CreateClusters. A list.
#' @param which.cluster Integer should correspond to the order of clusters
#' displayed.
#' @param fit.cluster How well should the wikipathways fit into the this
#' cluster? Integer from 0-1. Default is 0.99.
#'
#' @return A dataframe of which pathways corresponded best with the chosen
#' dynamics seen in the selected cluster.
#' @export
#'
#' @usage ReturnCluster(ClusterData, which.cluster, fit.cluster)
#' @examples
#' library(Mfuzz)
#' w <- list(WP1 = c("1071", "11303", "11806", "11808", "11812", "11813",
#' "11814", "11816", "13122", "13350", "15357", "15450",
#' "16816",  "16835", "16956", "16971", "17777", "18830",
#' "20652", "20778"),
#' WP10 = c("11651", "12702", "14784", "16186", "16198", "16199",
#' "16367", "16451", "16453", "18708", "19247", "20416",
#' "20846", "20848", "20850", "20851", "26395", "26396",
#' "26413", "26417",
#' "269523", "384783", "54721", "81601"),
#' WP103 = c("110196", "13121",  "13360", "14137", "15357", "16987",
#' "17855", "18194", "192156", "20775", "208715", "235293",
#' "319554", "66234", "68603"),
#' WP108 = c("107585",  "107869", "109079", "109815", "114679", "12916",
#' "13370", "13371",  "14080", "14281", "14775", "14776",
#' "14778", "16476", "18024",  "18033", "18986", "19697",
#' "19946", "20226", "20341", "20342",  "20363", "20364",
#' "20683", "20687", "20768", "211006", "214580",
#' "223776", "232223", "26462", "27361", "28042", "280621",
#' "433003", "50493", "50880", "619883", "625249", "625281",
#' "65967", "664969",  "69227", "71787", "71984", "72657",
#' "74777", "75420", "75512",  "80795", "9403"),
#' WP113 = c("12159", "12387", "12393", "12399"))
#' e <- list(D1 = c("1071", "11303", "11806", "11808", "223776", "232223",
#' "26462", "27361", "17855", "18194", "192156" ),
#' D2 = c("19946", "20226", "20341", "20342",  "20363", "20364",
#' "20683", "20687", "20768", "211006", "214580", "16971",
#' "17777", "18830"),
#' D3 = c("1071", "11303", "11806", "11808", "11812", "11813",
#' "11814", "11816", "13122", "13350", "20652", "20778"))
#' WikiMatrix(e_list = e, wp_list = w) -> wmat
#' TurnPercent(wikiMatrix = wmat, rowInt = 4) -> pmat
#' CreateClusters(method = 'c', Percent_matrix = pmat,
#' no.clusters = 3, Variance = 0)
#' ClusterCheck(Clusters = Clusters, W = FALSE)
#' Quickfuzz(Mfuzzdata = Mfuzzdata, Clusters = Clusters, W = FALSE)
#' ReturnCluster(ClusterData = ClusterData, which.cluster = 1,
#' fit.cluster = 0.5)
ReturnCluster <- function(ClusterData, which.cluster, fit.cluster = 0.99){
        as.data.frame(ClusterData) -> X
        X[which(X[,which.cluster] > fit.cluster),] -> singlecluster
return(singlecluster)
}
