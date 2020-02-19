#' @title miRDB_data
#' @description Produces a list of mRNA-miR predicted interactions from miRDB
#' for the species of interest from miRDB. Requires the tidyverse package.
#' @param miRDB Most recent predicted targets data file. Download from here
#'  http://mirdb.org/download.html.
#' @param species Species of interest e.g. hsa or mmu.
#' @param orgDB organism db file specific for species being looked into.
#' @return A dataframe of mRNA-miR interactons, and other information
#' associated with these interactions.
#' @export
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @import tidyverse
#' @usage miRDB_data(miRDB, species = '', orgDB)
#' @examples
#' library(org.Mm.eg.db)
#' library(tidyverse)
#' library(clusterProfiler)
#'data.frame(row.names = c("5151919", "5161682", "5151921", "5151922"),
#'"V1" = c("mmu-let-7a-1-3p", "mmu-let-7c-2-3p", "mmu-let-7k",
#'"mmu-miR-100-3p"),
#'"V2" = c("NM_144958", "NM_178648", "NM_133355", "NM_172405"),
#'"V3" = c(53.89190, 53.99680, 83.25980,64.69390)) -> miRDB_data
#' miRDB_data(miRDB = miRDB_data, species = 'mmu',
#' orgDB = org.Mm.eg.db) -> miRDB_mmu
miRDB_data <- function(miRDB, species, orgDB){
if (missing(miRDB)) stop('Input miRDB file.');
if (missing(species)) stop('e.g hsa or mmu.');
if (missing(orgDB)) stop('Input a speices library.');
miRDB_data <- miR <- NULL
names(miRDB) <- c('miR', 'Target', 'Score')
miRDB %>%
filter(str_detect(miR, species)) -> miRDB_s
miRDB_mRNA <- bitr(miRDB_s$Target, fromType = 'REFSEQ', toType = 'SYMBOL',
OrgDb = orgDB)
mirDB_merged <- merge(miRDB_s, miRDB_mRNA, by.x = 'Target', by.y = 'REFSEQ',
all = TRUE)
miRDB_df <- data.frame(miRDB_Interactions = paste(mirDB_merged$miR,
':', mirDB_merged$SYMBOL,sep = ''),
miRDB_miR = mirDB_merged$miR,
miRDB_mRNA = mirDB_merged$SYMBOL)
}
