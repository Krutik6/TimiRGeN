#' @title GenesLists
#' @description Will produce a list of nested dataframes. For combined
#' analysis method = c, and for separated analysis method = s. Combined
#' analysis colnames should be as follows 'timepoint.resultstype'.
#' GenesList will make new dataframes separated by 'timepoint.'.
#' Separated colnames should be as follows 'genetype_timepoint.resulttype'.
#' GenesList will make separate lists for 'genetype_'. Each of these lists
#' will have dataframes which have been made by separating by 'timepoint.'.
#' Make sure to follow colname nomenclature carefully.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param genetic_data If 'c', this should be a dataframe with miR and mRNA
#' information together.
#' @param timeString If 'c', this should be a commmon string representing
#' 'timepoint' e.g. for H.1, H.10, H.20, timeString = 'H'.
#' @param miR_data If 's', a dataframe of microRNA data. Rownames are genes
#' and colnames are like so; genetype_timepoint.resulttype.
#' timepoint.resulttype should be the same in mRNA and miR data.
#' @param mRNA_data If 's', a dataframe of mRNA data. Rownames are genes
#' and colnames are like so; genetype_timepoint.resulttype.
#' timepoint.resulttype should be the same in mRNA and miR data.
#' @return A list of dataframes separated by features in the column names which
#' can be stored in the metadata area of an MAE object.
#' @export
#' @importFrom stringr str_extract
#' @usage GenesList(method, genetic_data, timeString, miR_data, mRNA_data)
#' @examples
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' 
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' 
#' AddPrefix(gene_df = MAE@ExperimentList$miR,
#' prefixString = "miR" ) -> MAE@ExperimentList$miR_p
#' AddPrefix(gene_df = MAE@ExperimentList$mRNA, 
#' prefixString = "mRNA") -> MAE@ExperimentList$mRNA_p
#' 
#' GenesList(method = "s", miR_data = MAE@ExperimentList$miR_p,
#' mRNA_data = MAE@ExperimentList$mRNA_p) -> MAE@metadata$genelist
#' 
#' CombineGenes(miR_data = MAE@ExperimentList$miR, mRNA_data = 
#' MAE@ExperimentList$mRNA) -> MAE@ExperimentList$genetic_data
#' 
#' GenesList(method = 'c', genetic_data = MAE@ExperimentList$genetic_data,
#' timeString = 'D') -> MAE@metadata$genelist
GenesList <- function(method, genetic_data, timeString, miR_data, mRNA_data){
if (method == 'c') {
if (missing(genetic_data)) stop('Input combined miR and mRNA data.
Colnames structure should be
timepoint.resulttype.');
if (missing(timeString)) stop('Input timepoint. E.g. if colnames were D1,
D2, D3; then timeString should be D.');
colnames(genetic_data) <- gsub(x = colnames(genetic_data), pattern = timeString,
replacement = 'TP')
X <- lapply(split.default(genetic_data, sub("(TP\\d+).*", "\\1",
names(genetic_data))), as.list)
names(X) <- gsub(names(X), pattern = 'TP', replacement = timeString) 
L1 <- lapply(X, data.frame, stringsAsFactors = FALSE)
L2 <- lapply(L1, function(DF) {rownames(DF) <- rownames(genetic_data); DF})
L3 <- L2[gtools::mixedsort(names(L2))]
X <- lapply(L3, function(x) {
colnames(x) <- sub(x = colnames(x), pattern = 'TP',
replacement = timeString)
x}
)
return(X)
} else if (method == 's') {
if (missing(miR_data)) stop('Input miR data.
Colnames structure should be
genetype_timepoint.resulttype.');
if (missing(mRNA_data)) stop('Input mRNA data. Colnames structure should
be genetype_timepoint.resulttype.');
miR_data <- as.data.frame(miR_data)
mRNA_data <- as.data.frame(mRNA_data)
genedata <- list(miR_data = miR_data, mRNA_data = mRNA_data)
genes.split <- lapply(genedata, function(df)
{
Unigenes <- unique(str_extract(names(df),"\\S+\\."))
List <- lapply(Unigenes,function(name){return(df[,grep(name,names(df),
fixed=TRUE)])})
names(List) <- Unigenes
return(List)
})
} else{ stop('Please insert method c for combined analysis or s for
seperate analysis')}
}
