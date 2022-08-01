# The following experiments compare different algorithms for approximate inference

# Load packages

library("SGS")
library("usethis")

# packages for benchmark visualization

library("ggplot2")
library("RColorBrewer")
library("ggpubr")
library("cowplot")
library("gridExtra")
library("matrixStats")

###############################
## Test Different Dimensions ##
###############################

# Choose network parameters

N_var1 <- 50
N_Obs1 <- as.integer(N_var1/2) # assume half of the variables are observed
N_var2 <- 100
N_Obs2 <- as.integer(N_var2/2) # assume half of the variables are observed
N_var3 <- 150
N_Obs3 <- as.integer(N_var3/2) # assume half of the variables are observed
N_var4 <- 200
N_Obs4 <- as.integer(N_var4/2) # assume half of the variables are observed

# Choose sampling parameters

samplingMethods <- c("GS","SGS", "LBP")
N_rep=10
N_samples=100
N_nets=100


# Perform benchmark with multiple nets

resultTable1 <- benchmarkMultipleNets(N_var1, N_Obs1, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")
resultTable2 <- benchmarkMultipleNets(N_var2, N_Obs2, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")
resultTable3 <- benchmarkMultipleNets(N_var3, N_Obs3, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")
resultTable4 <- benchmarkMultipleNets(N_var4, N_Obs4, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")

# Visualize the results of the benchmark

SGS:::ProcessResultTableNRMSE(resultTable1)
SGS:::ProcessResultTableNRMSE(resultTable2)
SGS:::ProcessResultTableNRMSE(resultTable3)
SGS:::ProcessResultTableNRMSE(resultTable4)

# ProcessResultTable(resultTable1)
# ProcessResultTable(resultTable2)
# ProcessResultTable(resultTable3)


# store benchmark results
# saveRDS(list("resultTable1"=resultTable1,"resultTable2"=resultTable2,"resultTable3"=resultTable3,"resultTable4"=resultTable4), file = "benchmark_res/benchmarkDims.rds")
# 
# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkDims.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("50","100","150","200"), 
             c("Number of Nodes n=50","Number of Nodes n=100","Number of Nodes n=150","Number of Nodes n=200"), "Dims")


#########################
## Test Different Nets ##
#########################

# Choose network parameters

N_var5 <- 100
N_Obs5 <- as.integer(N_var5/2) # assume half of the variables are observed
N_var6 <- 100
N_Obs6 <- as.integer(N_var6/2) # assume half of the variables are observed
N_var7 <- 100
N_Obs7 <- as.integer(N_var7/2) # assume half of the variables are observed
N_var8 <- 100
N_Obs8 <- as.integer(N_var8/2) # assume half of the variables are observed

# Choose sampling parameters

samplingMethods <- c("GS","SGS", "LBP")
N_rep=5
N_samples=100
N_nets=100


# Perform benchmark with multiple nets

resultTable1 <- benchmarkMultipleNets(N_var5, N_Obs5, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")
resultTable2 <- benchmarkMultipleNets(N_var6, N_Obs6, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="interEr")
resultTable3 <- benchmarkMultipleNets(N_var7, N_Obs7, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="barabasi")
resultTable4 <- benchmarkMultipleNets(N_var8, N_Obs8, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="watts")


# store benchmark results
# saveRDS(list("resultTable1"=resultTable1,"resultTable2"=resultTable2,"resultTable3"=resultTable3,"resultTable4"=resultTable4), file = "benchmark_res/benchmarkNets.rds")
# 
# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkNets.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

# Visualize the results of the benchmark
#
# ProcessResultTableNRMSE(resultTable5)
# ProcessResultTableNRMSE(resultTable6)
# ProcessResultTableNRMSE(resultTable7)
# ProcessResultTableNRMSE(resultTable8)

# final plots 

SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("E-R","E-R Isl.","B-A","W-S"), 
             c("Erdős–Rényi graphs","Erdős–Rényi island graphs","Barabási-Albert graph","Watts-Strogatz graph"), "Nets")



#############################
## Test Different Sparcity ##
#############################

# Choose network parameters

N_var5 <- 50
N_Obs5 <- as.integer(N_var5/2) # assume half of the variables are observed
N_var6 <- 50
N_Obs6 <- as.integer(N_var6/2) # assume half of the variables are observed
N_var7 <- 50
N_Obs7 <- as.integer(N_var7/2) # assume half of the variables are observed
N_var8 <- 50
N_Obs8 <- as.integer(N_var8/2) # assume half of the variables are observed

# Choose sampling parameters

samplingMethods <- c("GS","SGS", "LBP")
N_rep=5
N_samples=100
N_nets=50


# Perform benchmark with multiple nets

resultTable1 <- benchmarkMultipleNets(N_var5, N_Obs5, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er",N_neighbours=1.4)
resultTable2 <- benchmarkMultipleNets(N_var6, N_Obs6, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er",N_neighbours=1.9)
resultTable3 <- benchmarkMultipleNets(N_var7, N_Obs7, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er",N_neighbours=2.35)
resultTable4 <- benchmarkMultipleNets(N_var8, N_Obs8, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er",N_neighbours=2.76)

# Visualize the results of the benchmark

ProcessResultTableNRMSE(resultTable1)
ProcessResultTableNRMSE(resultTable2)
ProcessResultTableNRMSE(resultTable3)
ProcessResultTableNRMSE(resultTable4)

# store benchmark results
# saveRDS(list("resultTable1"=resultTable1,"resultTable2"=resultTable2,"resultTable3"=resultTable3,"resultTable4"=resultTable4), file = "benchmark_res/benchmarkSparcity.rds")

# load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkSparcity.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("S=2","S=3","S=4","S=5"), 
             c("Markov Blanket size S=2","Markov Blanket size S=3","Markov Blanket size S=4","Markov Blanket size S=5"), "Sparcity")



#####################################
## Different Fractions of Evidence ##
#####################################

# Choose network parameters

N_var1 <- 50
N_Obs1 <- as.integer(N_var1/5) # assume half of the variables are observed
N_var2 <- 50
N_Obs2 <- as.integer(N_var2*2/5) # assume half of the variables are observed
N_var3 <- 50
N_Obs3 <- as.integer(N_var3*3/5) # assume half of the variables are observed
N_var4 <- 50
N_Obs4 <- as.integer(N_var4*4/5) # assume half of the variables are observed

# Choose sampling parameters

samplingMethods <- c("GS","SGS", "LBP")
N_rep=5
N_samples=100
N_nets=100


# Perform benchmark with multiple nets

resultTable1 <- benchmarkMultipleNets(N_var1, N_Obs1, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")
resultTable2 <- benchmarkMultipleNets(N_var2, N_Obs2, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")
resultTable3 <- benchmarkMultipleNets(N_var3, N_Obs3, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")
resultTable4 <- benchmarkMultipleNets(N_var4, N_Obs4, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er")

# Visualize the results of the benchmark

ProcessResultTableNRMSE(resultTable1)
ProcessResultTableNRMSE(resultTable2)
ProcessResultTableNRMSE(resultTable3)
ProcessResultTableNRMSE(resultTable4)

# store benchmark results
# saveRDS(list("resultTable1"=resultTable1,"resultTable2"=resultTable2,"resultTable3"=resultTable3,"resultTable4"=resultTable4), file = "benchmark_res/benchmarkFractions.rds")
# 
# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkFractions.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("f=0.2","f=0.4","f=0.6","f=0.8"), 
             c("Evidence Fraction f=0.2","Evidence Fraction f=0.4","Evidence Fraction f=0.6","Evidence Fraction f=0.8"), "Fractions")


#########################
## Multiple Categories ##
#########################

# Choose network parameters

N_var1 <- 50
N_Obs1 <- as.integer(N_var1/2) # assume half of the variables are observed
N_var2 <- 50
N_Obs2 <- as.integer(N_var2/2) # assume half of the variables are observed
N_var3 <- 50
N_Obs3 <- as.integer(N_var3/2) # assume half of the variables are observed
N_var4 <- 50
N_Obs4 <- as.integer(N_var4/2) # assume half of the variables are observed

# Choose sampling parameters

samplingMethods <- c("GS","SGS", "LBP")
N_rep=5
N_samples=100
N_nets=50


# Perform benchmark with multiple nets

resultTable1 <- benchmarkMultipleNets(N_var1, N_Obs1, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er", nodeDim=2)
resultTable2 <- benchmarkMultipleNets(N_var2, N_Obs2, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er", nodeDim=4)
resultTable3 <- benchmarkMultipleNets(N_var3, N_Obs3, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er", nodeDim=6)
resultTable4 <- benchmarkMultipleNets(N_var4, N_Obs4, N_nets, N_rep, N_samples, samplingMethods, DAGmethod="er", nodeDim=8)

# Visualize the results of the benchmark

ProcessResultTableNRMSE(resultTable1)
ProcessResultTableNRMSE(resultTable2)
ProcessResultTableNRMSE(resultTable3)
ProcessResultTableNRMSE(resultTable4)

# ProcessResultTable(resultTable1)
# ProcessResultTable(resultTable2)
# ProcessResultTable(resultTable3)


# store benchmark results 
# saveRDS(list("resultTable1"=resultTable1,"resultTable2"=resultTable2,"resultTable3"=resultTable3,"resultTable4"=resultTable4), file = "benchmark_res/benchmarkCategories.rds")
# 
# load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkCategories.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("C=2","C=4","C=6", "C=8"), 
             c("Categories C=2","Categories C=4","Categories C=6","Categories C=8"), "Categories")








