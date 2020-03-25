#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(ggplot2)
#load data
sigwiki <- readRDS("EnrichWiki.rds")
MAE <- MultiAssayExperiment()
metadata(MAE)[["sigwiki"]] <- sigwiki
#test Quickbar
p1 <- Quickbar(X = metadata(MAE)[[1]][[1]],Y = metadata(MAE)[[1]][1])
#test plot
p2 <- ggplot(head(sigwiki$D7_wikipathways[which(
                                 sigwiki$D7_wikipathways$p.adjust < 0.05),],
                                 n = 15), 
                                 aes(x=reorder(Description, -p.adjust),
                                 y=Count, 
                                 fill=-p.adjust)) +
    geom_bar(stat = "identity", width = 0.5) +
    
    scale_fill_continuous(type = "gradient") +
    
    labs(y = "Associated genes", x = "wikipathways", fill = "p.adjust") +
    
    theme(axis.text=element_text(size=18)) +
    
    ggtitle(names(sigwiki[4])) +
    
    theme(plot.title = element_text(2, face = "bold", hjust = 0.5),
                                    legend.key.size = unit(2, "line")) +
    
    theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
    
    coord_flip()
#visual check
test_that("should be 4 wikipathways", {
    expect_equal(length(rownames(sigwiki$D7_wikipathways@result)),4)})
#outputs are the same


