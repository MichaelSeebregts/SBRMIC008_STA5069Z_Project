library(plotly)
library(ggplot2)
library(reticulate)
use_condaenv("py310", required = TRUE)

library(tensorflow)

library(keras)
library(lle)

# GMM's 
library(mclust)
library(dbscan)

set.seed(56)

pixelsMatrixRaw = readRDS("F:/.University/5th Year/STA5069Z/SBRMIC008_STA5069Z_Project/pixelsMatrix")
patientInfoRaw = read.csv("F:/.University/5th Year/STA5069Z/SBRMIC008_STA5069Z_Project/PatientInfo.csv")

obsNotNA = which(!is.na(patientInfoRaw$cdr))

patientInfo = patientInfoRaw[obsNotNA, ]

imagesWithNACDR = sprintf("%04d", obsNotNA)

# Removing 2 cdr and grouping it with 1 cdr
patientInfo$cdr[patientInfo$cdr == 2] = 1


front = paste0("image_",imagesWithNACDR, "_front")
top = paste0("image_",imagesWithNACDR, "_top")
side = paste0("image_",imagesWithNACDR, "_side")

pixelsMatrix = list()

for (i in 1:length(obsNotNA))
{
  tempFront = front[i]
  tempTop = top[i]
  tempSide = side[i]
  
  pixelsMatrix[[tempFront]] = pixelsMatrixRaw[[tempFront]]
  pixelsMatrix[[tempSide]] = pixelsMatrixRaw[[tempSide]]
  pixelsMatrix[[tempTop]] = pixelsMatrixRaw[[tempTop]]
  
}

flatPixels = lapply(pixelsMatrix, c)

flatFront = matrix(flatPixels[["image_0001_front"]], nrow = 1)
flatTop = matrix(flatPixels[["image_0001_top"]], nrow = 1)
flatSide = matrix(flatPixels[["image_0001_side"]], nrow = 1)


for (j in 2:length(obsNotNA))
{
  
  flatFront = rbind(flatFront, t(flatPixels[[front[j]]]))
  flatTop = rbind(flatTop, t(flatPixels[[top[j]]]))
  flatSide = rbind(flatSide, t(flatPixels[[side[j]]]))
}

frontDF = data.frame(flatFront, cdr = as.factor(patientInfo$cdr))
topDF = data.frame(flatTop, cdr = as.factor(patientInfo$cdr))
sideDF = data.frame(flatSide, cdr = as.factor(patientInfo$cdr))

frontPCA = prcomp(flatFront, center = TRUE, scale = FALSE)
# frontPCADF = data.frame("PC1" = frontPCA$x[, 1], "PC2" = frontPCA$x[, 2], "cdr" = frontDF$cdr)
# ggplot(frontPCADF, aes(x = PC1, y = PC2, color = cdr)) + geom_point()

topPCA = prcomp(flatTop, center = TRUE, scale = FALSE)
# topPCADF = data.frame("PC1" = topPCA$x[, 1], "PC2" = topPCA$x[, 2], "cdr" = topDF$cdr)
# ggplot(topPCADF, aes(x = PC1, y = PC2, color = cdr)) + geom_point()

sidePCA = prcomp(flatSide, center = TRUE, scale = FALSE)
# sidePCADF = data.frame("PC1" = sidePCA$x[, 1], "PC2" = sidePCA$x[, 2], "cdr" = sideDF$cdr)
# ggplot(sidePCADF, aes(x = PC1, y = PC2, color = cdr)) + geom_point()

cumulativeVarianceFront = cumsum(summary(frontPCA)$importance[2, ])
plot(cumulativeVarianceFront)

cumulativeVarianceTop = cumsum(summary(topPCA)$importance[2, ])
plot(cumulativeVarianceTop)

cumulativeVarianceSide = cumsum(summary(sidePCA)$importance[2, ])
plot(cumulativeVarianceSide)

kmeansFront = kmeans(frontPCA$x[, 1:100], centers = 3)
kmeansSide = kmeans(topPCA$x[, 1:100], centers = 3)
kmeansTop = kmeans(sidePCA$x[, 1:100], centers = 3)

PCADFFront = data.frame("clusters" = kmeansFront$cluster, "cdr" = patientInfo$cdr)
PCADFSide = data.frame("clusters" = kmeansSide$cluster, "cdr" = patientInfo$cdr)
PCADFTop = data.frame("clusters" = kmeansTop$cluster, "cdr" = patientInfo$cdr)

table(PCADFFront)
table(PCADFSide)
table(PCADFTop)

kmeansFront2 = kmeans(frontPCA$x[, 1:100], centers = 2)
kmeansSide2 = kmeans(topPCA$x[, 1:100], centers = 2)
kmeansTop2 = kmeans(sidePCA$x[, 1:100], centers = 2)

PCADFFront2 = data.frame("clusters" = kmeansFront2$cluster, "cdr" = patientInfo$cdr)
PCADFSide2 = data.frame("clusters" = kmeansSide2$cluster, "cdr" = patientInfo$cdr)
PCADFTop2 = data.frame("clusters" = kmeansTop2$cluster, "cdr" = patientInfo$cdr)

table(PCADFFront2)
table(PCADFSide2)
table(PCADFTop2)

optFront = optics(frontPCA$x[,1:100], minPts = 50)
plot(optFront)

dbscanFront = extractDBSCAN(optFront, eps_cl = 20)
PCADBScanDF = data.frame("clusters" = dbscanFront$cluster, "cdr" = patientInfo$cdr)
table(PCADBScanDF)
# ggplot(frontPCADF, aes(x = PC1, y = PC2, color = as.factor(kmeansFront$cluster))) + 
#   geom_text(aes(label = cdr))
# ggplot(topPCADF, aes(x = PC1, y = PC2, color = as.factor(kmeansTop$cluster))) + 
#   geom_text(aes(label = cdr))
# ggplot(sidePCADF, aes(x = PC1, y = PC2, color = as.factor(kmeansSide$cluster))) + 
#   geom_text(aes(label = cdr))

# LLE 

lleResult = lle(X = flatFront, m = 70, k = 30)
lleComp = lleResult$Y

head(lleComp)

kmeansLLE = kmeans(lleComp, centers = 3)
dfLLE = data.frame("Cluster" = kmeansLLE$cluster, "cdr" = patientInfo$cdr)
table(dfLLE)

kmeansLLE2 = kmeans(lleComp, centers = 2)
dfLLE2 = data.frame("Cluster" = kmeansLLE2$cluster, "cdr" = patientInfo$cdr)
table(dfLLE2)

optLLE = optics(lleComp, minPts = 50)
plot(optLLE)

dbscanLLE = extractDBSCAN(optLLE, eps_cl = 10)
LLEDBScanDF = data.frame("clusters" = dbscanLLE$cluster, "cdr" = patientInfo$cdr)
table(LLEDBScanDF)

# Auto Encoder for Front Scans
modelFront = keras_model_sequential() %>% 
  layer_dense(units = 2056, activation = "linear", input_shape = ncol(flatFront)) %>% 
  layer_dense(units = 1024, activation = "sigmoid") %>% 
  layer_dense(units = 128, activation = "sigmoid", name = "bottleneck") %>% 
  layer_dense(units = 1024, activation = "sigmoid") %>% 
  layer_dense(units = 2056, activation = "linear") %>% 
  layer_dense(units = ncol(flatFront))

modelFront %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelFront %>% fit(x = flatFront, 
              y = flatFront, 
              epochs = 25, 
              batch_size = 128, 
              validation_split = 0.2, 
              verbose = 0)

modelFront2 = keras_model_sequential() %>% 
  layer_dense(units = 2056, activation = "relu", input_shape = ncol(flatFront)) %>% 
  layer_dense(units = 1024, activation = "relu") %>% 
  layer_dense(units = 128, activation = "relu", name = "bottleneck") %>% 
  layer_dense(units = 1024, activation = "relu") %>% 
  layer_dense(units = 2056, activation = "relu") %>% 
  layer_dense(units = ncol(flatFront))

modelFront2 %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelFront2 %>% fit(x = flatFront, 
                   y = flatFront, 
                   epochs = 25, 
                   batch_size = 128, 
                   validation_split = 0.2, 
                   verbose = 0)

intermediateLayerFront = keras_model(inputs = modelFront$input, 
                                     outputs = get_layer(modelFront, "bottleneck")$output)

intermediateOutputFront = predict(intermediateLayerFront, flatFront)

autoEncFrontKmeans = kmeans(intermediateOutputFront, centers = 3)

AutoEncFrontDF = data.frame("clusters" = autoEncFrontKmeans$cluster, "cdr" = patientInfo$cdr)
table(AutoEncFrontDF)

autoEncFrontKmeans2 = kmeans(intermediateOutputFront, centers = 2)

AutoEncFrontDF2 = data.frame("clusters" = autoEncFrontKmeans2$cluster, "cdr" = patientInfo$cdr)
table(AutoEncFrontDF2)

optFrontAE = optics(intermediateOutputFront, minPts = 50)
plot(optFrontAE)

dbscanFrontAE = extractDBSCAN(optFrontAE, eps_cl = 10)
AEFrontDBScanDF = data.frame("clusters" = dbscanFrontAE$cluster, "cdr" = patientInfo$cdr)
table(AEFrontDBScanDF)

# Auto Encoder for Top Scans
modelTop = keras_model_sequential() %>% 
  layer_dense(units = 2056, activation = "linear", input_shape = ncol(flatTop)) %>% 
  layer_dense(units = 1024, activation = "sigmoid") %>% 
  layer_dense(units = 128, activation = "sigmoid", name = "bottleneck") %>% 
  layer_dense(units = 1024, activation = "sigmoid") %>% 
  layer_dense(units = 2056, activation = "linear") %>% 
  layer_dense(units = ncol(flatTop))

modelTop %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelTop %>% fit(x = flatTop, 
                   y = flatTop, 
                   epochs = 25, 
                   batch_size = 128, 
                   validation_split = 0.2, 
                   verbose = 0)

modelTop2 = keras_model_sequential() %>% 
  layer_dense(units = 2056, activation = "relu", input_shape = ncol(flatTop)) %>% 
  layer_dense(units = 1024, activation = "relu") %>% 
  layer_dense(units = 128, activation = "relu", name = "bottleneck") %>% 
  layer_dense(units = 1024, activation = "relu") %>% 
  layer_dense(units = 2056, activation = "relu") %>% 
  layer_dense(units = ncol(flatTop))

modelTop2 %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelTop2 %>% fit(x = flatTop, 
                 y = flatTop, 
                 epochs = 25, 
                 batch_size = 128, 
                 validation_split = 0.2, 
                 verbose = 0)

intermediateLayerTop = keras_model(inputs = modelTop$input, 
                                     outputs = get_layer(modelTop, "bottleneck")$output)

intermediateOutputTop = predict(intermediateLayerTop, flatTop)

autoEncTopKmeans = kmeans(intermediateOutputTop, centers = 3)

AutoEncTopDF = data.frame("clusters" = autoEncTopKmeans$cluster, "cdr" = patientInfo$cdr)
table(AutoEncTopDF)

autoEncTopKmeans2 = kmeans(intermediateOutputTop, centers = 2)

AutoEncTopDF2 = data.frame("clusters" = autoEncTopKmeans2$cluster, "cdr" = patientInfo$cdr)
table(AutoEncTopDF2)

optTopAE = optics(intermediateOutputTop, minPts = 50)
plot(optTopAE)

dbscanTopAE = extractDBSCAN(optTopAE, eps_cl = 10)
AETopDBScanDF = data.frame("clusters" = dbscanTopAE$cluster, "cdr" = patientInfo$cdr)
table(AETopDBScanDF)


# Auto Encoder for Side Scans
modelSide = keras_model_sequential() %>% 
  layer_dense(units = 2056, activation = "linear", input_shape = ncol(flatSide)) %>% 
  layer_dense(units = 1024, activation = "sigmoid") %>% 
  layer_dense(units = 128, activation = "sigmoid", name = "bottleneck") %>% 
  layer_dense(units = 1024, activation = "sigmoid") %>% 
  layer_dense(units = 2056, activation = "linear") %>% 
  layer_dense(units = ncol(flatSide))

modelSide %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelSide %>% fit(x = flatSide, 
                   y = flatSide, 
                   epochs = 25, 
                   batch_size = 128, 
                   validation_split = 0.2, 
                   verbose = 0)


modelSide2 = keras_model_sequential() %>% 
  layer_dense(units = 2056, activation = "relu", input_shape = ncol(flatSide)) %>% 
  layer_dense(units = 1024, activation = "relu") %>% 
  layer_dense(units = 128, activation = "relu", name = "bottleneck") %>% 
  layer_dense(units = 1024, activation = "relu") %>% 
  layer_dense(units = 2056, activation = "relu") %>% 
  layer_dense(units = ncol(flatSide))

modelSide2 %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelSide2 %>% fit(x = flatSide, 
                  y = flatSide, 
                  epochs = 25, 
                  batch_size = 128, 
                  validation_split = 0.2, 
                  verbose = 0)

intermediateLayerSide = keras_model(inputs = modelSide$input, 
                                     outputs = get_layer(modelSide, "bottleneck")$output)

intermediateOutputSide = predict(intermediateLayerSide, flatSide)

autoEncKSidemeans = kmeans(intermediateOutputSide, centers = 3)

AutoEncSideDF = data.frame("clusters" = autoEncKSidemeans$cluster, "cdr" = patientInfo$cdr)
table(AutoEncSideDF)

autoEncKSidemeans2 = kmeans(intermediateOutputSide, centers = 2)

AutoEncSideDF2 = data.frame("clusters" = autoEncKSidemeans2$cluster, "cdr" = patientInfo$cdr)
table(AutoEncSideDF2)

optSideAE = optics(intermediateOutputSide, minPts = 50)
plot(optSideAE)

dbscanSideAE = extractDBSCAN(optSideAE, eps_cl = 10)
AESideDBScanDF = data.frame("clusters" = dbscanSideAE$cluster, "cdr" = patientInfo$cdr)
table(AESideDBScanDF)

# Convolution for Front Scan

numImages = nrow(flatFront)
imgFrontHeight = dim(pixelsMatrix$image_0001_front)[1]
imgFrontWidth = dim(pixelsMatrix$image_0001_front)[2]

flatFrontReshaped = array(flatFront, dim = c(numImages, imgFrontHeight, imgFrontWidth, 1))
flatFrontReshaped = aperm(flatFrontReshaped, c(1, 3, 2, 4), resize = FALSE)

modelConvFront = keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "linear", padding = "same", 
                input_shape = c(imgFrontHeight, imgFrontWidth, 1)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same") %>%
  
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu", padding = "same") %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same") %>%
  
  layer_conv_2d(filters = 4, kernel_size = c(2, 2), activation = "relu", padding = "same") %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same", name = "bottleneck") %>%
  
  # Decoder
  layer_conv_2d(filters = 4, kernel_size = c(2, 2), activation = "relu", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "linear", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 1, kernel_size = c(3, 3), activation = "linear", padding = "same")

  
  
modelConvFront %>% compile(
  loss = "mse", 
  optimizer = "adam"
  )

modelConvFront %>% fit(x = flatFrontReshaped,
                   y = flatFrontReshaped,
                   epochs = 50,
                   batch_size = 16,
                   validation_split = 0.2,
                   verbose = 0)

mseFConv = evaluate(modelConvFront, flatFrontReshaped, flatFrontReshaped)
mseFConv

intermediateLayerFrontConv = keras_model(inputs = modelConvFront$input,
                                     outputs = get_layer(modelConvFront, "bottleneck")$output)

intermediateOutputFrontConv = predict(intermediateLayerFrontConv, flatFrontReshaped)


bottleneckFlat = array_reshape(intermediateOutputFrontConv, c(dim(intermediateOutputFrontConv)[1], -1))

kmeansFrontConv = kmeans(bottleneckFlat, centers = 3)

convFrontDFClust = data.frame("cluster" = kmeansFrontConv$cluster, "cdr" = patientInfo$cdr)

table(convFrontDFClust)

kmeansFrontConv2 = kmeans(bottleneckFlat, centers = 2)

convFrontDFClust2 = data.frame("cluster" = kmeansFrontConv2$cluster, "cdr" = patientInfo$cdr)

table(convFrontDFClust2)

optFrontConv = optics(bottleneckFlat, minPts = 50)
plot(optFrontConv)

dbscanFrontConv = extractDBSCAN(optFrontConv, eps_cl = 10)
ConvFrontDBScanDF = data.frame("clusters" = dbscanFrontConv$cluster, "cdr" = patientInfo$cdr)
table(ConvFrontDBScanDF)

# Convolution for Top Scan

numImages = nrow(flatTop)
imgTopHeight = dim(pixelsMatrix$image_0001_top)[1]
imgTopWidth = dim(pixelsMatrix$image_0001_top)[2]

flatTopReshaped = array(flatTop, dim = c(numImages, imgTopHeight, imgTopWidth, 1))
# flatTopReshaped = aperm(flatTopReshaped, c(1, 3, 2, 4), resize = FALSE)

modelConvTop = keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "linear", padding = "same", 
                input_shape = c(imgTopHeight, imgTopWidth, 1)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same") %>%
  
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu", padding = "same") %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same") %>%
  
  layer_conv_2d(filters = 4, kernel_size = c(2, 2), activation = "relu", padding = "same") %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same", name = "bottleneck") %>%
  
  # Decoder
  layer_conv_2d(filters = 4, kernel_size = c(2, 2), activation = "relu", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "linear", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 1, kernel_size = c(3, 3), activation = "linear", padding = "same")



modelConvTop %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelConvTop %>% fit(x = flatTopReshaped,
                  y = flatTopReshaped,
                  epochs = 50,
                  batch_size = 16,
                  validation_split = 0.2,
                  verbose = 0)

mseTopConv = evaluate(modelConvTop, flatTopReshaped, flatTopReshaped)
mseTopConv

intermediateLayerTopConv = keras_model(inputs = modelConvTop$input,
                                         outputs = get_layer(modelConvTop, "bottleneck")$output)

intermediateOutputTopConv = predict(intermediateLayerTopConv, flatTopReshaped)


bottleneckTopFlat = array_reshape(intermediateOutputTopConv, c(dim(intermediateOutputTopConv)[1], -1))

kmeansTopConv = kmeans(bottleneckTopFlat, centers = 3)

convTopDFClust = data.frame("cluster" = kmeansTopConv$cluster, "cdr" = patientInfo$cdr)

table(convTopDFClust)

kmeansTopConv2 = kmeans(bottleneckTopFlat, centers = 2)

convTopDFClust2 = data.frame("cluster" = kmeansTopConv2$cluster, "cdr" = patientInfo$cdr)

table(convTopDFClust2)

optTopConv = optics(bottleneckTopFlat, minPts = 50)
plot(optTopConv)

dbscanTopConv = extractDBSCAN(optTopConv, eps_cl = 10)
ConvTopDBScanDF = data.frame("clusters" = dbscanTopConv$cluster, "cdr" = patientInfo$cdr)
table(ConvTopDBScanDF)

# Convolution for Side Scan

numImages = nrow(flatSide)
imgSideHeight = dim(pixelsMatrix$image_0001_side)[1]
imgSideWidth = dim(pixelsMatrix$image_0001_side)[2]

sideMatrix = matrix(flatSide, nrow = numImages, byrow = FALSE)

flatSideReshaped = array(sideMatrix, dim = c(numImages, imgSideHeight, imgSideWidth, 1))
flatSideReshaped = aperm(flatSideReshaped, c(1, 3, 2, 4))

# identical(t(flatSideReshaped[1,,,]), pixelsMatrix$image_0001_side)

modelConvSide = keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "linear", padding = "same", 
                input_shape = c(imgSideWidth, imgSideHeight, 1)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same") %>%
  
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu", padding = "same") %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same") %>%
  
  layer_conv_2d(filters = 4, kernel_size = c(2, 2), activation = "relu", padding = "same") %>%
  layer_max_pooling_2d(pool_size = c(2, 2), padding = "same", name = "bottleneck") %>%
  
  # Decoder
  layer_conv_2d(filters = 4, kernel_size = c(2, 2), activation = "relu", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "linear", padding = "same") %>%
  layer_upsampling_2d(size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 1, kernel_size = c(3, 3), activation = "linear", padding = "same")



modelConvSide %>% compile(
  loss = "mse", 
  optimizer = "adam"
)

modelConvSide %>% fit(x = flatSideReshaped,
                  y = flatSideReshaped,
                  epochs = 50,
                  batch_size = 16,
                  validation_split = 0.2,
                  verbose = 0)

mseSideConv = evaluate(modelConvSide, flatSideReshaped, flatSideReshaped)
mseSideConv

intermediateLayerSideConv = keras_model(inputs = modelConvSide$input,
                                         outputs = get_layer(modelConvSide, "bottleneck")$output)

intermediateOutputSideConv = predict(intermediateLayerSideConv, flatSideReshaped)


bottleneckSide = array_reshape(intermediateOutputSideConv, c(dim(intermediateOutputSideConv)[1], -1))

kmeansSideConv = kmeans(bottleneckSide, centers = 3)

convSideDFClust = data.frame("cluster" = kmeansSideConv$cluster, "cdr" = patientInfo$cdr)


table(convSideDFClust)

kmeansSideConv2 = kmeans(bottleneckSide, centers = 2)

convSideDFClust2 = data.frame("cluster" = kmeansSideConv2$cluster, "cdr" = patientInfo$cdr)


table(convSideDFClust)

optSideConv = optics(bottleneckSide, minPts = 50)
plot(optSideConv)

dbscanSideConv = extractDBSCAN(optSideConv, eps_cl = 10)
ConvSideDBScanDF = data.frame("clusters" = dbscanSideConv$cluster, "cdr" = patientInfo$cdr)
table(ConvSideDBScanDF)

# Full Data 

combinedScans = cbind(flatFront, flatTop, flatSide)


