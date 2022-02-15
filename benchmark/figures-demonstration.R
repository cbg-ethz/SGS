# The following experiments compare different algorithms for approximate inference

# Load packages

library("usethis")
library("devtools")
library("pcalg")
library("ggplot2")
library("RColorBrewer")
library("ggpubr")

install_github("PhysFritz/SubGroupSeparation")
# load_all()

# packages for benchmark visualization

library("cowplot")
library("gridExtra")
library("matrixStats")

###############################################################################
## Test Different Dimensions ##################################################
###############################################################################

# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkDims.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("50","100","150","200"), 
             c("Number of Nodes N=50","Number of Nodes N=100","Number of Nodes N=150","Number of Nodes N=200"), "Dims2", width = 6, height = 4)

###############################################################################
## Test Different Nets ########################################################
###############################################################################

# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkNets.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

# final plots 

makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("E-R","E-R Isl.","B-A","W-S"), 
             c("Erdoes-Renyi graph","Erdoes-Renyi Island graph","Barabási-Albert graph","Watts-Strogatz graph"), "Nets2", width = 6, height = 4)

###############################################################################
## Test Different Sparcity ####################################################
###############################################################################

# load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkSparcity.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("S=2","S=3","S=4","S=5"), 
             c("Markov Blanket size S=2","Markov Blanket size S=3","Markov Blanket size S=4","Markov Blanket size S=5"), "Sparcity2", width = 6, height = 4)

###############################################################################
## Different Fractions of Evidence ############################################
###############################################################################

# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkFractions.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("f=0.2","f=0.4","f=0.6","f=0.8"), 
             c("Evidence Fraction f=0.2","Evidence Fraction f=0.4","Evidence Fraction f=0.6","Evidence Fraction f=0.8"), "Fractions2", width = 6, height = 4)

###############################################################################
## Multiple Categories ########################################################
###############################################################################

# load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkCategories.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("C=2","C=4","C=6", "C=8"), 
             c("Categories C=2","Categories C=4","Categories C=6","Categories C=8"), "Categories2", width = 6, height = 4)




