setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(igraph)
#load data
filt_df <- structure(list(avecor = c(-0.929199786400515, -0.729228501795928,
-0.431983639087243, -0.55088842103792, -0.978422379116014,
-0.627856061946295, -0.998864281242427, -0.877362389481312, -0.990125035397228,
-0.771338310408749), miR = structure(c(9L, 5L, 6L, 2L, 8L, 4L, 3L, 4L,
7L, 1L), .Label = c("hsa-miR-107", "hsa-miR-193a-3p", "hsa-miR-28-5p",
"hsa-miR-331-3p", "hsa-miR-362-3p", "hsa-miR-362-5p", "hsa-miR-429",
"hsa-miR-590-5p", "hsa-miR-630" ), class = "factor"),
mRNA = structure(c(1L, 2L, 2L, 3L, 3L, 4L, 5L, 5L,
5L, 6L), .Label = c("IGF1R", "PRKCA", "TESK2", "THBS1",  "TLN2", "VAV3"),
class = "factor"), miR_Entrez = structure(c(8L,  7L, 7L, 4L, 6L, 3L, 5L,
3L, 1L, 2L), .Label = c("ENSG00000198976",  "ENSG00000198997",
"ENSG00000199172", "ENSG00000207614", "ENSG00000207651",
"ENSG00000207741", "ENSG00000208015", "ENSG00000283798"),
class = "factor"), mRNA_Entrez = structure(c(4L, 5L, 5L, 1L, 1L, 3L, 6L,
6L, 6L, 2L), .Label = c("ENSG00000070759", "ENSG00000134215",
"ENSG00000137801", "ENSG00000140443", "ENSG00000154229","ENSG00000171914"),
class = "factor"), TargetScan = c(0L,  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L), miRDB = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
Predicted_Interactions = c(1L,  2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L),
miRTarBase = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), Pred_Fun = c(2L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)), class = "data.frame",
row.names = c("hsa-miR-630:IGF1R", "hsa-miR-362-3p:PRKCA",
"hsa-miR-362-5p:PRKCA", "hsa-miR-193a-3p:TESK2", "hsa-miR-590-5p:TESK2",
"hsa-miR-331-3p:THBS1", "hsa-miR-28-5p:TLN2", "hsa-miR-331-3p:TLN2",
"hsa-miR-429:TLN2", "hsa-miR-107:VAV3" ))
#test function
MakeNet(filt_df = filt_df) -> net
#check 1
test_that("net is an igraph object", {
expect_equal(class(net), "igraph")
expect_equal(length(net), 10)
})
#internal checks
rbind(filt_df, filt_df) -> df
df$id <- "id"
genes <- data.frame(genes = c(as.character(filt_df$miR),
as.character(filt_df$mRNA)))
#check 2
test_that("20 values in genes", {
expect_equal(length(rownames(genes)), 20)
})
#continue
as.integer(factor(genes[,1])) -> genes$id
paste("s", genes$id, sep = "") -> df$id
rownames(df)<- NULL
max(as.integer(rownames(df))/2) -> halfway
#check 3
test_that("10 values in halfway", {
expect_equal(halfway, 10)
})
#continue
nodes <- data.frame(id = df$id,
genes = c(as.character(df$miR[seq_len(halfway)]),
as.character(df$mRNA[seq_len(halfway)])))
nodes$genetype <- "mRNA"
max(as.integer(rownames(nodes))/2) -> number_miRs
for (i in seq_len(number_miRs)) {
nodes$genetype[i] <- "miR"
}
#check 4
test_that("aspects of nodes", {
expect_equal(length(rownames(nodes)), 20)
expect_equal(length(colnames(nodes)), 3)
})
#continue
links <- data.frame(from = nodes$id[seq_len(halfway)],
to = nodes$id[-c(seq_len(halfway))],
Databases = filt_df$Pred_Fun,
Correlation = filt_df$avecor,
type = "hyperlink")
#check 5
test_that("aspects of links", {
expect_equal(length(rownames(links)), 10)
expect_equal(length(colnames(links)), 5)
})
#continue
nodes[! duplicated(nodes$genes),] -> nodes1
net1 <- graph_from_data_frame(d=links, vertices=nodes1, directed=TRUE)
plot(net1)
dev.off()
#check 6
test_that("aspects of net", {
expect_equal(length(net), length(net1))
})
#save data
saveRDS(net, "MakeNet/net.rds")
