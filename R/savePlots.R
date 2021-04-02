#' @title savePlots
#' @description Saves all plots from enrichWiki into the current working
#' directory.
#' @param largeList A large list containing GSEA results. This should be
#' stored as metadata within the MAE used in the enrichWiki function.
#' @param maxInt Integer, number of samples in data set.
#' @param fileType Type of file for images to be exported as: "png", "tiff",
#' "svg" or "jpeg".
#' @param width = Width of plots in inches. Default is 22 inches.
#' @param height = Heightof plots in inches. Default is 10 inches.
#' @return Saves plots in working directory. Each sample (e.g. time point) will
#' have a separate plot.
#' @export
#' @importFrom grDevices colorRampPalette dev.off jpeg png svg tiff x11
#' @usage savePlots(largeList, maxInt, fileType = '', width, height)
savePlots <- function(largeList, maxInt, fileType, width=22, height=13){

    if (missing(largeList)) stop('largeList is missing. Add large list of nested dataframe. This is the output from the enrichWiki function and should be stored as metadata of the MAE used in the enrichWiki function.')

    if (missing(maxInt)) stop('maxInt is missing. Add number of samples your data has. Should be an integer')

  # create empty list
  plot_list <- list()

  # create plot object
  for (i in seq_len(maxInt)) {
    p <- quickBar(X = largeList[[i]]@result, Y = names(largeList[i]))
    plot_list[[i]] = p

  }

  # create file names
  for (i in seq_len(maxInt)) {
    file_name <- paste(names(largeList[i]), ".", fileType, sep="")

    # options for plotting

    if (fileType == "tiff") {
      tiff(file_name, width = width, height = height, units = 'in', res = 100)
      print(plot_list[[i]])
      dev.off()

    } else if (fileType == "png") {
      png(file_name, width = width, height = height, units = 'in', res = 100)
      print(plot_list[[i]])
      dev.off()

    } else if (fileType == "svg") {
      svg(file_name, width = width, height = height)
      print(plot_list[[i]])
      dev.off()

    } else if (fileType == "jpeg") {
      jpeg(file_name, width = width, height = height, units = 'in', res = 100)
      print(plot_list[[i]])
      dev.off()

    } else {print("Input relevant file type for export: tiff, png, svg or jpeg")}
  }
}
