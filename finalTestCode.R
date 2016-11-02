# Final Test Code:
library(DGSA)

# setwd("~/Google Drive/Functional_Research/DGSA")
# source("./R/main.R")

INPUT         <- read.csv("./data/Input.csv")
OUTPUT.damage <- read.csv("./data/Output_damage.csv")

comps  <- prcomp(OUTPUT.damage)
scores <- t(t(comps$x[,1:2])/comps$sdev[1:2])
clustering = kmeans(scores, 2, 50)$cluster
plot(scores)
points(scores[clustering == 1,], col="red", pch = 19)
points(scores[clustering == 2,], col="blue", pch=19)

# EDA for DGSA:
plotCDFS(clustering, INPUT, .code = "all*")
plotCDFS(clustering, INPUT, .code = "beta")
plotCDFS(clustering, INPUT, .code = "lambda")
plotCDFS(clustering, INPUT, .code = "beta|lambda", .nBins = 2)
plotCDFS(clustering, INPUT, .code = "lambda|beta", .nBins = 2)


# run DGSA:
system.time(myDGSA <- dgsa(clustering, INPUT, .normalize = TRUE, .interactions = TRUE, .nBoot = 100))
plotMatrixDGSA(myDGSA)
plotMatrixDGSA(myDGSA, .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single')
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single', .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single', tl.srt = 65)
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single', .hypothesis = FALSE, tl.srt = 65)

plotMatrixDGSA(myDGSA, .method = "circle", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "number", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "square", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "shade", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "pie", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "pie", type = "lower", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "pie", type = "upper", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "circle", diag=TRUE, .hypothesis = TRUE, tl.srt = 45,
               insig = "pch", pch = "+", pch.col = "black", pch.cex = 1.5)

plotParetoDGSA(myDGSA)
plotParetoDGSA(myDGSA, .interaction = "beta")
plotParetoDGSA(myDGSA, .interaction = "lambda")



# Trees are much more robust in a sense that you are averaging over all splits and all trees.
# Therefore sensitivities are always consistent. In dgsa they are not!


# Shale Dataset:
library(DGSA)
load("~/Google Drive/Functional_Research/TreeS_workingDir/TreeS_workingDir/RealDataSetScaled.Rdata")

CovCompletion <- Reservoir1_allCovariates[,colnames(covariatesCompletion)]
CovCompletion <- cbind(CovCompletion,Reservoir1_allCovariates$GeolX_Rel, Reservoir1_allCovariates$GeolY_Rel )
colnames(CovCompletion) <- c(colnames(covariatesCompletion), "GeolX_Rel", "GeolY_Rel")

comps  <- prcomp(t(completeSet))
scores <- t(t(comps$x[,1:2])/comps$sdev[1:2])
clustering = kmeans(scores, 2, 50)$cluster
plot(scores)
points(scores[clustering == 1,], col="red", pch = 19)
points(scores[clustering == 2,], col="blue", pch=19)

# EDA for DGSA:

# colnames(CovCompletion)
# [1] "Geol3DSpacing"                 "GeolAPIGrav"                   "ProdVert200"
# [4] "ProdVert300"                   "ProdVertVin200"                "ProdVertVin300"
# [7] "ComplStagesPumped"             "ComplStageInterval"            "ComplFluid"
# [10] "ComplProppant"                 "ComplAvgRate"                  "ComplAvgPress"
# [13] "DrillWellTort"                 "CMP_LATERAL_LENGTH"            "CMP_STIMULATED_LATERAL_LENGTH"
# [16] "CMP_STAGES_STIMULATED"         "CMP_AVERAGE_INTERVAL"          "CMP_AVG_BREAK_PRESSURE"
# [19] "CMP_TOTAL_FLUID_STIM"          "CMP_CLEAN_FLUID_TOTAL"         "CMP_TOTAL_PROPPANT_USED"
# [22] "CMP_BAKE_TIME.days."           "GeolX_Rel"                     "GeolY_Rel"
#
CovCompletion <- CovCompletion[,-which(colnames(CovCompletion) %in% c("ProdVertVin200", "ProdVertVin300", "CMP_AVG_BREAK_PRESSURE"))]
plotCDFS(clustering, CovCompletion, .code = "all*")
plotCDFS(clustering, CovCompletion, .code = "GeolAPIGrav")
plotCDFS(clustering, CovCompletion, .code = "CMP_STIMULATED_LATERAL_LENGTH")
plotCDFS(clustering, CovCompletion, .code = "CMP_STIMULATED_LATERAL_LENGTH|GeolX_Rel", .nBins = 2)
plotCDFS(clustering, CovCompletion, .code = "ProdVertVin200|CMP_BAKE_TIME.days.", .nBins = 2)
plotCDFS(clustering, CovCompletion, .code = "CMP_STIMULATED_LATERAL_LENGTH|ProdVertVin200", .nBins = 2)
plotCDFS(clustering, CovCompletion, .code = "ProdVert200|CMP_STIMULATED_LATERAL_LENGTH", .nBins = 2)
plotCDFS(clustering, CovCompletion, .code = "ProdVertVin200")
plotCDFS(clustering, CovCompletion, .code = "ProdVert200")
plotCDFS(clustering, CovCompletion, .code = "GeolY_Rel")
plotCDFS(clustering, CovCompletion, .code = "GeolAPIGrav")
plotCDFS(clustering, CovCompletion, .code = "GeolAPIGrav|GeolY_Rel")
plotCDFS(clustering, CovCompletion, .code = "GeolY_Rel|GeolAPIGrav")

# run DGSA:
myDGSA <- DGSA::dgsa(clustering, CovCompletion, .interactions = TRUE, .nBoot = 100)
plotMatrixDGSA(myDGSA)
plotMatrixDGSA(myDGSA, .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single')
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single', .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single', tl.srt = 65)
plotMatrixDGSA(myDGSA, order = 'hclust', hclust.method='single', .hypothesis = FALSE, tl.srt = 65)

plotMatrixDGSA(myDGSA, .method = "circle", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "number", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "square", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "shade", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "pie", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "pie", type = "lower", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "pie", type = "upper", .hypothesis = FALSE)
plotMatrixDGSA(myDGSA, .method = "circle", order = 'hclust', hclust.method='centroid',
               diag=TRUE, .hypothesis = TRUE, tl.srt = 90,
               insig = "pch", pch = "*", pch.col = "blue", pch.cex = 1.5)


plotCDFS(clustering, CovCompletion, .code = "CMP_STIMULATED_LATERAL_LENGTH|GeolAPIGrav", .nBins = 3)
plotCDFS(clustering, CovCompletion, .code = "GeolAPIGrav|CMP_STIMULATED_LATERAL_LENGTH", .nBins = 3)

plotCDFS(clustering, CovCompletion, .code = "GeolAPIGrav|GeolY_Rel", .nBins = 3)
plotCDFS(clustering, CovCompletion, .code = "GeolY_Rel|GeolAPIGrav", .nBins = 3)


plotCDFS(clustering, CovCompletion, .code = "GeolY_Rel", .nBins = 2)
plotCDFS(clustering, CovCompletion, .code = "GeolAPIGrav", .nBins = 2)

plot(CovCompletion$GeolY_Rel, CovCompletion$GeolAPIGrav,
     main = paste("Correlation Coef:" ,round(cor(CovCompletion$GeolY_Rel, CovCompletion$GeolAPIGrav),2)))

plotParetoDGSA(myDGSA)
plotParetoDGSA(myDGSA, .interaction = "beta")
plotParetoDGSA(myDGSA, .interaction = "lambda")


### Boostrapping by considering clustering as well.
INPUT         <- read.csv("./data/Input.csv")
OUTPUT.damage <- read.csv("./data/Output_damage.csv")

comps  <- prcomp(OUTPUT.damage)
scores <- t(t(comps$x[,1:2])/comps$sdev[1:2])
clustering = kmeans(scores, 2, 50)$cluster
baseDGSA <- dgsa(clustering, INPUT, .normalize = FALSE, .interactions = TRUE, .nBoot = 100, .progress = TRUE)

library(doParallel)
library(foreach)

cl <- makeCluster(8)
registerDoParallel(cl)

bootList <- foreach(c.seed = 1:500, .packages = c('DGSA')) %dopar% {

  set.seed(c.seed)
  sampleInd <- sample(nrow(INPUT), nrow(INPUT), replace=TRUE)


  dgsa(clustering[sampleInd], INPUT[sampleInd,], .normalize = FALSE,.interactions = TRUE, .nBoot = 100, .progress=FALSE)

}
stopCluster(cl)

bootExtractCL1 <- sapply(bootList, function(x) x$sensitivityMatrix[1,,], simplify = 'array')
bootExtractCL2 <- sapply(bootList, function(x) x$sensitivityMatrix[2,,], simplify = 'array')

cl1Normalizer  <- apply(bootExtractCL1, c(1,2), function(x){
  if(any(is.na(x))){
    return(NaN)
  } else {
    quantile(x, 0.95)
  }
}
)

cl2Normalizer  <- apply(bootExtractCL2, c(1,2), function(x){
  if(any(is.na(x))){
    return(NaN)
  } else {
    quantile(x, 0.95)
  }
}
)

normalizerMatrix      <- baseDGSA$sensitivityMatrix
normalizerMatrix[1,,] <- cl1Normalizer
normalizerMatrix[2,,] <- cl2Normalizer

normalizedSensitivity <- baseDGSA
normalizedSensitivity$sensitivityMatrix <- normalizedSensitivity$sensitivityMatrix / normalizerMatrix
originalDGSA <- dgsa(clustering, INPUT, .normalize = TRUE, .interactions = TRUE, .nBoot = 500, .progress = TRUE)
class(normalizedSensitivity) <- class(originalDGSA)

plotMatrixDGSA(originalDGSA, .hypothesis = FALSE)
plotMatrixDGSA(normalizedSensitivity, .hypothesis = FALSE)



