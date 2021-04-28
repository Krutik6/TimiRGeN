#' @title clusterList
#' @description clusterList will transform clusters created by createClusters2
#' into lists based on which genes associate most to each cluster. Genes which
#' associate with a cluster are determined by the fitCluster parameter in the
#' function.
#' @param MAE MultiAssayExperiment which will store the results from
#' createClusters.
#' @param clusterData A dataframe which contains cluster-pathway fit scores
#' and is stored as an assay within the MAE used in the createClusters2 function.
#' @param fitCluster Integer from 0-1. How well should genes fit into a
#' cluster? Default is 0.5.
#' @param miR_IDs miR_ensembl or miR_entrez. Use a getIDsMir function to
#' acquire this. This will be stored as an assay in the MAE used in a getIdsMir
#' function.
#' @param mRNA_IDs mRNA_ensembl or mRNA_entrez. Use a getIDsMrna function to
#' acquire this. This will be stored as an assay in the MAE used in a getIdsMrna
#' function.
#' @return A list containing the genes which fit to each cluster.
#' @export
#' @usage clusterList(MAE, clusterData, fitCluster, miR_IDs, mRNA_IDs)
#' @examples
#' library(org.Hs.eg.db)
#'
#' data(long_data)
#'
#' miRNA <- long_data[c(1:105),]
#'
#' mRNA <- long_data[-c(1:105),]
#'
#' MAE <- startObject(miR = miRNA, mRNA = mRNA)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Hs.eg.db, 'hsa')
#'
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), mirror = 'useast', 'hsapiens',
#'                   orgDB = org.Hs.eg.db)
#'
#' MAE <- combineGenes(MAE, miR_data = assay(MAE, 1),
#'                     mRNA_data = assay(MAE, 2))
#'
#' MAE <- createClusters2(MAE = MAE, genetic_data = assay(MAE, 9),
#'                        noClusters =2)
#'
#' MAE <- clusterList(MAE = MAE, clusterData = assay(MAE, 11), fitCluster = 0.5,
#'                    miR_IDs = assay(MAE, 3),
#'                    mRNA_IDs = assay(MAE, 7))
clusterList <- function(MAE, clusterData, fitCluster = 0.5, miR_IDs, mRNA_IDs) {

  if (missing(MAE)) stop('MAE is missing. Add MAE. Please use createClusters2 first.')

  if (missing(clusterData)) stop('clusterData is missing. Add dataframe which has cluster-pathway fit scores. Please use the createClusters2 function first. clusterData should be stored as an assay within the MAE used in the createClusters2 function.')

  if (missing(miR_IDs)) stop('miR_IDs is missing. Add dataframe of miR gene IDs. Please use getIdsMir first. Output of a getIdsMir function should be stored as assays within the MAE used in the getIdsMir function.')

  if (missing(mRNA_IDs)) stop('mRNA_IDs is missing. Add dataframe of mRNA gene IDs. Please use getIds first. Output of getIds function should be stored as an assay within the MAE used in the getIdsMrna function.')

  metadata <- `metadata<-` <- NULL

  Clist <- list()

  for (i in seq_along(1:ncol(clusterData))) {

    Clist[[i]] <- clusterData[which(clusterData[,i] > fitCluster),]

    names(Clist)[[i]] <- paste0("Cluster ",names(clusterData)[[i]])

  }

  MAE2 <- MultiAssayExperiment()

  MAE2 <- TimiRGeN::addIds(MAE = MAE2, method = "c",
                           filtered_genelist = Clist,
                           miR_IDs = miR_IDs,
                           mRNA_IDs = mRNA_IDs)

  MAE2 <- TimiRGeN::eNames(MAE = MAE2, method = "c",
                           gene_IDs = metadata(MAE2)[[1]])

  ID_list <- metadata(MAE2)[[2]]

  metadata(MAE)[["ID_list"]] <- ID_list

  return(MAE)
}
