#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(igraph)
#load data
net <- readRDS("net.rds")
Quicknet(metadata(net)[[1]])
dev.off()
#visual check, legends a bit wonky
test_that("net should be 10 long", {
  expect_equal(length(metadata(net)[[1]]), 10)
})
