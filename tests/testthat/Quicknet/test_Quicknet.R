setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(igraph)
#load data
readRDS("MakeNet/net.rds") -> net
Quicknet(net = net)
#visual check, legends a bit wonky
test_that("net should be 10 long", {
  expect_equal(length(net), 10)
})
