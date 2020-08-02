#' @title savePlots
#' @description Saves all plots from enrichWiki into the current directory.
#' @param largeList A large list containing GSEA results. This should be
#' contained as metadata within the MAE used in the enrichWiki function.
#' @param maxInt Integer, number of samples in data set.
#' @param quickType quickDot or quickBar. This will be the plot type.
#' @param fileType Type of file for images to be exported as: "png", "tiff",
#' "svg" or "jpg."
#' @return saves plots in current directory.
#' @export
#' @importFrom grDevices colorRampPalette dev.off jpeg png svg tiff x11
#' @usage savePlots(largeList, maxInt, quickType, fileType = '')
#' @examples
#' \donttest{
#' library(org.Mm.eg.db)
#'
#' MAE <- MultiAssayExperiment()
#'
#' metadata(MAE)[["e_list"]] <- e_list
#'
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
#'
#' MAE <- enrichWiki(MAE = MAE, method = 'c', ID_list = metadata(MAE)[[1]],
#'                  orgDB = org.Mm.eg.db, path_gene = assay(MAE, 1),
#'                  path_name = assay(MAE, 2), ID = "ENTREZID",
#'                  universe = assay(MAE, 1)[[2]])
#'
#' savePlots(largeList = metadata(MAE)[[2]], maxInt = 5,
#'          quickType = quickDot, fileType = "jpg")
#
#'}
savePlots <- function(largeList, maxInt, quickType, fileType){

    if (missing(largeList)) stop('Add large list of nested dataframe. This
                                 is the output from the enrichWiki function and
                                 should be stored as metadata of the MAE used
                                 in the enrichWiki function.')

    if (missing(maxInt)) stop('Add number of samples your data has.
                              Should be an integer')

    if (missing(quickType)) stop('Add type of plot. quickBar for a bar plot
                                 or quickDot for a dot plot.')

    if (missing(fileType)) stop('Input the type of file to be saved either
                                 "png", "tiff", "svg" or "jpg".')

    # create empty list
    plot_list <- list()

    # create plot object
    for (i in seq_len(maxInt)) {
        p <- quickType(X = largeList[[i]]@result, Y = largeList[i])
        plot_list[[i]] = p

    }

    # create file names
    for (i in seq_len(maxInt)) {
        file_name <- paste(names(largeList[i]), ".", fileType, sep="")

    # options for plotting

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
