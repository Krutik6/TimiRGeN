#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(igraph)
#load data
readRDS("net.rds") -> net
Quicknet(net = net)
dev.off()
#visual check, legends a bit wonky
test_that("net should be 10 long", {
  expect_equal(length(net), 10)
})
