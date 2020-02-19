#' @title SavePlots
#' @description Saves all possible plots in the current directory.
#' @param largeList A list of wikipathways associated to entrezIDs for each
#' sample.
#' @param maxInt Number of samples in dataset/ lists in sigwiki
#' @param quickType Quickplot ot Quickbar
#' @param fileType Type of file for images to be exported as: png, tiff,
#' svg or jpg
#' @return saves tiff files of all plots in current dir
#' @export
#' @importFrom grDevices colorRampPalette dev.off jpeg png svg tiff x11
#' @usage SavePlots(largeList, maxInt, quickType, fileType = '')
#' @examples
#' library(biomaRt)
#' library(ggplot2)
#' library(org.Mm.eg.db)
#' library(clusterProfiler)
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' 
#' e_list -> MAE@metadata$elist
#' dloadGMT(MAE, speciesInitial = "Mm") -> MAE
#' 
#' MAE@metadata$sigwiki <- EnrichWiki(method = "c",
#' e_list = MAE@metadata$elist,
#' orgDB = org.Mm.eg.db, 
#' path_gene = MAE@ExperimentList$path_gene, 
#' path_name = MAE@ExperimentList$path_name, 
#' ID = "ENTREZID", 
#' universe = MAE@ExperimentList$path_gene$gene)
#' 
#' SavePlots(largeList = MAE@metadata$sigwiki, maxInt = 5, 
#' quickType = Quickdot, fileType = "jpg") -> P
SavePlots <- function(largeList, maxInt, quickType, fileType){
if (missing(largeList)) stop('Input large list of nested dataframe.');
if (missing(maxInt)) stop('Input number of samples your data has.');
if (missing(quickType)) stop('Input type of plot.');
if (missing(fileType)) stop('Input the type of file to be saved either
png, tiff, svg or jpg')
plot_list <- list()
for (i in seq_len(maxInt)) {
p <- quickType(X = largeList[[i]]@result, Y = largeList[i])
plot_list[[i]] = p
}
for (i in seq_len(maxInt)) {
file_name <- paste(names(largeList[i]), ".", fileType, sep="")
if (fileType == "tiff") {
tiff(file_name, width = 985, height = 847)
print(plot_list[[i]])
dev.off()
} else if (fileType == "png") {
png(file_name, width = 985, height = 847)
print(plot_list[[i]])
dev.off()
} else if (fileType == "svg") {
svg(file_name)
print(plot_list[[i]])
dev.off()
} else if (fileType == "jpg") {
jpeg(file_name, width = 985, height = 847)
print(plot_list[[i]])
dev.off()
} else ("Input relevant file type for export: tiff, png, svg or jpg")
}
}

