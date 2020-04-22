#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(igraph)
#load data
filt_df <- readRDS("filt_df.rds")
MAE <- MultiAssayExperiment()
#test function
net <- makeNet(MAE, filt_df = filt_df)
#check 1
test_that("net is an igraph object", {
    expect_equal(class(metadata(net)[[1]]), "igraph")
    expect_equal(length(metadata(net)[[1]]), 10)
})
#internal checks
df <- rbind(filt_df, filt_df)
df$id <- "id"
genes <- data.frame(genes = c(as.character(filt_df$miR),
                    as.character(filt_df$mRNA)))
#check 2
test_that("20 values in genes", {
    expect_equal(length(rownames(genes)), 10)
})
#continue
genes$id <- as.integer(factor(genes[,1]))
df$id <- paste("s", genes$id, sep = "")
rownames(df)<- NULL
halfway <- max(as.integer(rownames(df))/2)
#check 3
test_that("10 values in halfway", {
    expect_equal(halfway, 5)
})
#continue
nodes <- data.frame(id = df$id,
                    genes = c(as.character(df$miR[seq_len(halfway)]),
                              as.character(df$mRNA[seq_len(halfway)])))
nodes$genetype <- "mRNA"
number_miRs <- max(as.integer(rownames(nodes))/2)

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
nodes1 <- nodes[! duplicated(nodes$genes),]
net1 <- graph_from_data_frame(d=links, vertices=nodes1, directed=TRUE)
plot(net1)
dev.off()
#check 6
test_that("aspects of net", {
    expect_equal(length(metadata(net)[[1]]), length(net1))
})
metadata(MAE)[["net"]] <- net
#save data
saveRDS(net, "net.rds")
