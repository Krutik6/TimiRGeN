#' @title dloadmiRTarBase
#' @description Downloads most recent version (7.0) of functional targets from
#'miRTarBase database.
#' @return Dataframe of functionally tested miR-mRNA interactions in different
#' species.
#' @export
#' @usage dloadmiRTarBase()
dloadmiRTarBase <- function(){
download.file(
'http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/miRTarBase_MTI.xlsx',
'miRTarBase.xlsx')
readxl::read_excel('miRTarBase.xlsx') -> miRTarBase
write.csv(miRTarBase, 'miRTarBase.csv')
read.csv('miRTarBase.csv', row.names = 1) -> miRTarBase
print("Downloaded published miRNA target interaction data version 7.0")
print("Check if a newer version is available?
http://mirtarbase.mbc.nctu.edu.tw/php/index.php." )
print("If so download that and run it through miRTarBase_data function.")
file.remove("miRTarBase.xlsx")
return(miRTarBase)
}
