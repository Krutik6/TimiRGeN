#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(igraph)
#load data
readRDS("filt_df.rds") -> filt_df
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
expect_equal(length(rownames(genes)), 10)
})
#continue
as.integer(factor(genes[,1])) -> genes$id
paste("s", genes$id, sep = "") -> df$id
rownames(df)<- NULL
max(as.integer(rownames(df))/2) -> halfway
#check 3
test_that("10 values in halfway", {
expect_equal(halfway, 5)
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
expect_equal(length(rownames(nodes)), 10)
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
expect_equal(length(rownames(links)), 5)
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
saveRDS(net, "net.rds")
