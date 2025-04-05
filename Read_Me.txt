Important note, not all the data could be uploaded to GitHub due to file size restrictions. The pixels matrix can be sent upon request. Raw data is available on the OASIS1 website but was unable to be uploaded due to file size restrictions. 

This project is investigating the use of multivariate statistical techniques in the classification of MRI scans into different levels of Alzheimers. The MRI scans were first dimension reduced and then clustering techniques were applied to the reduced dimensions to see whether it would provide accurate predictions. A Convolutional Neural Network was used as well and its results compared to that of the dimensions reduction and clustering techniques. 


Data Dictionary: 

pixelsMatrix - File containing all the pixels values for the top, front and side scans of each individual in a matrix. The values are stored as a list 

PatientInfo - CSV file containing all the additional information provided for an individual. 

obsNotNA - the position of each individual without a CDR rating 

imagesWithNACDR - finds the individuals MRI scans for individuals without a CDR rating

flatFront - Converting the matrix of pixel values into a vector for each individuals front scan who have a CDR rating

flatTop - Converting the matrix of pixel values into a vector for each individuals top scan who have a CDR rating

flatSide- Converting the matrix of pixel values into a vector for each individuals side scan who have a CDR rating

frontPCA/topPCA/sidePCA - Principal Components for the front, top and side scans

kmeansFront3LLE/kmeansFront3AE - K-Means clustering for LLE and AutoEncoders with 3 clusters for the front scan. Equivalent are made for top and side scans

kmeansFront3LLE/kmeansFront2AE - K-Means clustering for LLE and AutoEncoders with 2 clusters for the front scan. Equivalent are made for top and side scans

optLLEFront/optAEFront - using the optics package on the LLE and AutoEncoder for the front scan. Equivalent are made for top and side scans

DBScanLLEFront/ DBScanAEFront - Using DBScan for LLE optics and AutoEncoder optics for the front scan. Equivalent are maded for top and sided scans. 

modelFront/ model2Front - 2 different AutoEncoder models for front scans. Equivalent made for side and top scans. 

mseFrontAE/ mseFrontAE2 - MSE for modelFront and modelFront2. Equivalent made for side and top scans. 

flatFrontReshaped/ flatTopReshaped/ flatSideReshaped - Reshaping the vector of number of images, imageheight, image width and color scheme in order to feed it into the CNN. 

input_layerFront/ input_layerTop/ input_layerSide - shape of the input layer for the CNN

branch1Front/ branch2Front/ branch3Front - 3 branches which will be added together and fed into the nodes. Equivalent made for top and side scans. 

mergedFront/ mergedTop/ mergedSide - Merged 3 branches for front, top and side scans. 

outputFront/ outputTop/ outputSide - Using merged 3 branches and feeding into the 128, 64 and 3 nodes. 

multiBranchCNNFront/ multiBranchCNNTop/ multiBranchCNNSide - combining the entire CNN

mseFConv/ mseTConv/ mseSConv - MSE for front, side and top scans. 

predictionsClassFront/ predictionsClassTop/ predictionsClassSide - predictions for each of the inputs for front, top and side scans. 

misClassifiedFront/ misClassifiedTop/ misClassifiedSide - Checking which predictions are misclassified for front, top and side scans. 
