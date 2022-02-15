# reproducable R script of application``

# load BiDAG packages

library(BiDAG)
library(reshape2)
library(pcalg)
library(ggplot2)
library(Bestie)

# Load packages

library("usethis")
library("devtools")
library("RColorBrewer")
library("ggpubr")

library("parallel")

# load_all()
install_github("PhysFritz/SubGroupSeparation")

################################
# Application replication code # 
################################


# import the Korean kidney cancer data
mutation_array <- readRDS(file = "mutation_array_processed.rds")
patienet_data <- readRDS(file = "patient_array_processed.rds")


## Step 1: Learn the full Bayes nets of each cancer subtype ##
##############################################################

set.seed(40)

# set structure learning learning parameters
myBdepar <- list("chi" = 4, "edgepf" = 1)

#constructing score objects
edgepmat <- string2mat(colnames(kirp), interactions, mapping, type = "pf", pf = 2)
KIRPscore <- scoreparameters("bde", kirp, bdepar = myBdepar,
                               edgepmat = edgepmat)
KIRCscore <- scoreparameters("bde", kirc, bdepar = myBdepar,
                               edgepmat = edgepmat)
  
#Estimating MAP DAGs

#KIRP
itKirp <- iterativeMCMC(KIRPscore, mergetype = "skeleton", verbose = FALSE)
#KIRC
itKirc <- iterativeMCMC(KIRCscore, mergetype = "skeleton", verbose = FALSE)

# get the condtional probability tables (CPTs) of both graphs
kircDAGaugmented <- Bestie:::DAGparameters(itKirc$DAG, KIRCscore)
kirpDAGaugmented <- Bestie:::DAGparameters(itKirp$DAG, KIRPscore)

# convert Bayes net of Bestie to SGS format
kircBNSGS <- convertBNBestieToSGS(kircDAGaugmented)
kirpBNSGS <- convertBNBestieToSGS(kirpDAGaugmented)


## Step 2: Classification of Korean data via normalizing constant ##
####################################################################

# get mutual genes of both cohorts
mutualGenes <- intersect(colnames(kirc),colnames(mutation_array)) # 26 (of 70 and 184)
# get position of mutual genes in Bayes net and mutation array
mutualGenesInBN <- sapply(mutualGenes, function(xx) which(colnames(kirc) == xx))
mutualGenesInMA <- sapply(mutualGenes, function(xx) which(colnames(mutation_array) == xx))

n_patients <- length(patienet_data[,1])
normConstSample <- array(NA, dim=c(n_patients,4))

  
# for each patient, get the corresponding mutation value of mutual genes
patientObs <- list()
n_patients <- length(patienet_data[,1])
for (tempPatient in 1:n_patients){
  # create observation for each patient for Bayesian network inference
  tempObs <- list()
  tempObs$observed.vars <- mutualGenesInBN
  tempObs$observed.vals <- mutation_array[tempPatient,mutualGenesInMA]+1
  # in the observation of the binary varible, "+1" accounts for the fact
  # that SGSs binary scale is "1,2" and mutations are as "0,1"
  patientObs[[tempPatient]] <- tempObs
}

# get the normalizing constant for each patient for the two networks
normConstSample <- array(NA, dim=c(n_patients,4))
rownames(normConstSample) <- patienet_data[,1]
colnames(normConstSample) <- c("NC_KIRC","NC_KIRP","Assigned.Cancer.Type","Real.Cancer.Type")
# add the real cancer type as a reference
normConstSample[,4] <- patienet_data[,14]
for (tempPatient2 in 1:length(patientObs)){
  # show progress
  print(paste("Sample" ,tempPatient2, "of", length(patientObs)))
  NCkirc <- junction.tree.normConstSGS(kircBNSGS, patientObs[[tempPatient2]])
  NCkirp <- junction.tree.normConstSGS(kirpBNSGS, patientObs[[tempPatient2]])
  normConstSample[tempPatient2,1:2] <- c(NCkirc,NCkirp)
  if(NCkirc>NCkirp){
    normConstSample[tempPatient2,3] <- "ccRCC"
  }else{
    normConstSample[tempPatient2,3] <- "pRCC"
  }
}

# save
inferenceRes <- normConstSample


# now repeat step 1 and 2, but without the normalizing constant, removing all
# all missing varibles from the analysis

## Step 3: Learn the reduced Bayes nets of each cancer subtype ##
#################################################################

# select data with mutual genes only
kirpReduced <- kirp[,mutualGenesInBN]
kircReduced <- kirc[,mutualGenesInBN]

set.seed(17)
  
#constructing score objects 
edgepmat_mutual <- string2mat(colnames(kirpReduced), interactions, mapping, type = "pf", pf = 2)
KIRPscore_mutual <- scoreparameters("bde", kirpReduced, bdepar = myBdepar,
                                    edgepmat = edgepmat_mutual)
KIRCscore_mutual <- scoreparameters("bde", kircReduced, bdepar = myBdepar,
                                    edgepmat = edgepmat_mutual)

#Estimating MAP DAGs

#KIRP
itKirp_mutual <- iterativeMCMC(KIRPscore_mutual, mergetype = "skeleton", verbose = FALSE)
#KIRC
itKirc_mutual <- iterativeMCMC(KIRCscore_mutual, mergetype = "skeleton", verbose = FALSE)

# get the condtional probability tables (CPTs) of both graphs
kirpDAGaugmented_mutual <- Bestie:::DAGparameters(itKirp_mutual$DAG, KIRPscore_mutual)
kircDAGaugmented_mutual <- Bestie:::DAGparameters(itKirc_mutual$DAG, KIRCscore_mutual)

# convert Bayes net of Bestie to SGS format
kircBNSGS_mutual <- convertBNBestieToSGS(kircDAGaugmented_mutual)
kirpBNSGS_mutual <- convertBNBestieToSGS(kirpDAGaugmented_mutual)

# get mutual genes of both cohorts
# mutualGenes <- intersect(rownames(itKirp$DAG),colnames(mutation_array)) # 26 (of 70 and 184)


########################################################################
## Step 4: Classification of Korean data without normalizing constant ##
########################################################################
  
# get position of mutual genes in Bayes net and mutation array
mutualGenesInBN <- sapply(mutualGenes, function(xx) which(colnames(kirpReduced) == xx))
mutualGenesInMA <- sapply(mutualGenes, function(xx) which(colnames(mutation_array) == xx))

# for each patient, get the corresponding mutation value of mutual genes
patientObs <- list()
n_patients <- length(mutation_array[,1])
for (tempPatient in 1:n_patients){
  # create observation for each patient for Bayesian network inference
  tempObs <- list()
  tempObs$observed.vars <- mutualGenesInBN
  tempObs$observed.vals <- mutation_array[tempPatient,mutualGenesInMA]+1
  # in the observation of the binary varible, "+1" accounts for the fact
  # that SGSs binary scale is "1,2" and mutations are as "0,1"
  patientObs[[tempPatient]] <- tempObs
}

# get the normalizing constant for each patient for the two networks
normConstSample <- array(NA, dim=c(n_patients,4))
rownames(normConstSample) <- rownames(mutation_array)
colnames(normConstSample) <- c("NC_KIRC","NC_KIRP","Assigned.Cancer.Type","Real.Cancer.Type")
# add the real cancer type
normConstSample[,4] <- patienet_data[,14]

for (tempPatient2 in 1:length(mutation_array[,1])){
  NCkirc <- junction.tree.normConstSGS(kircBNSGS_mutual, patientObs[[tempPatient2]])
  NCkirp <- junction.tree.normConstSGS(kirpBNSGS_mutual, patientObs[[tempPatient2]])
  normConstSample[tempPatient2,1:2] <- c(NCkirc,NCkirp)
  if(NCkirc>NCkirp){
    normConstSample[tempPatient2,3] <- "ccRCC"
  }else{
    normConstSample[tempPatient2,3] <- "pRCC"
  }
}

# store the Bayes nets
inferenceRes_mutual <- normConstSample

# summarize results
results <- list("inferenceRes"=inferenceRes,"inferenceRes_mutual"=inferenceRes_mutual)

# saveRDS(results, "analysis_results.rds")
# results <- readRDS(results, "analysis_results")

###############################
## Step 5: Visualise results ##
###############################

# create barplot
# with normalizing constant
plot.barplotSamplesKPS(results$inferenceRes)
# without normalizing constant
plot.barplotSamplesKPS(results$inferenceRes_mutual)

# create barplot including the difference
plot.barplotSamplesKPS_difference(results$inferenceRes,results$inferenceRes_mutual)

# plot ROC curve
plot.ROC_curve_KPS(results$inferenceRes,results$inferenceRes_mutual)



