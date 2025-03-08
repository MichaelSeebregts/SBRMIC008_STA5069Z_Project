library(magick)
library(EBImage)
library(stringr)
library(dplyr)

# Install BiocManager to install EBImage
#install.packages("BiocManager", dependencies = TRUE)
# Use BiocManager to install EBImage
#BiocManager::install("EBImage")

# Following code was used to read in all the images, convert them into their pixels
# These matrices of pixels were normalised and saved to a list denoting if they are top, front, side
# This list was written to files in order to not have to read in the raw images everytime. 

# listOfFiles = data.frame(c())
# listOfPixels = list()
# 
# main = "F:/.University/5th Year/STA5069Z/SBRMIC008_STA5069Z_Project/Data/Oasis/"
# 
# file_names = list.files(main, full.names = TRUE)
# 
# numID = c()
# 
# for (i in 1:length(file_names))
# {
#   numID = c(numID, (str_extract(file_names[i], "(?<=OAS1_)\\d{4}")))
#   sub = paste(file_names[i],"/PROCESSED/MPRAGE/T88_111", sep = "")
# 
#   front = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n4_anon_111_t88_gfc_cor_110.gif", sep = "")
#   front2 = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n3_anon_111_t88_gfc_cor_110.gif", sep = "")
#   front3 = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n6_anon_111_t88_gfc_cor_110.gif", sep = "")
#   if (file.exists(front))
#   {
#     side = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n4_anon_111_t88_gfc_sag_95.gif", sep = "")
#     top = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n4_anon_111_t88_gfc_tra_90.gif", sep = "")
#     tempObs = c(front, side, top)
#     listOfFiles = rbind(listOfFiles, tempObs)
#   } else if (file.exists(front2))
#   {
#     side = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n3_anon_111_t88_gfc_sag_95.gif", sep = "")
#     top = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n3_anon_111_t88_gfc_tra_90.gif", sep = "")
#     tempObs = c(front2, side, top)
#     listOfFiles = rbind(listOfFiles, tempObs)
#   }else if (file.exists(front3))
#   {
#     side = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n6_anon_111_t88_gfc_sag_95.gif", sep = "")
#     top = paste(sub,"/OAS1_",numID[i],"_MR1_mpr_n6_anon_111_t88_gfc_tra_90.gif", sep = "")
#     tempObs = c(front3, side, top)
#     listOfFiles = rbind(listOfFiles, tempObs)
#   }
# 
# }
# image_list = list()
# 
# for (j in 1:length(file_names))
# {
#   imgFront = image_read(listOfFiles[j, 1])
#   imgFront = image_convert(imgFront, colorspace = "gray")
#   imgFront_dat = image_data(imgFront)
# 
#   imgSide = image_read(listOfFiles[j, 2])
#   imgSide = image_convert(imgSide, colorspace = "gray")
#   imgSide_dat = image_data(imgSide)
# 
#   imgTop = image_read(listOfFiles[j, 3])
#   imgTop = image_convert(imgTop, colorspace = "gray")
#   imgTop_dat = image_data(imgTop)
# 
#   heightFront = dim(imgFront_dat)[3]
#   widthFront = dim(imgFront_dat)[2]
# 
#   heightSide = dim(imgSide_dat)[3]
#   widthSide = dim(imgSide_dat)[2]
# 
#   heightTop = dim(imgTop_dat)[3]
#   widthTop = dim(imgTop_dat)[2]
# 
#   full_imgFront_dat = matrix(as.numeric(imgFront_dat[,,1]), ncol = 1, byrow = TRUE)
#   full_imgSide_dat = matrix(as.numeric(imgSide_dat[,,1]), ncol = 1, byrow = TRUE)
#   full_imgTop_dat = matrix(as.numeric(imgTop_dat[,,1]), ncol = 1, byrow = TRUE)
# 
#   for (q in 2:heightFront)
#   {
#     full_imgFront_dat = cbind(full_imgFront_dat, as.numeric(imgFront_dat[,,q]))
#   }
# 
#   for (z in 2:heightSide)
#   {
#     full_imgSide_dat = cbind(full_imgSide_dat, as.numeric(imgSide_dat[,,z]))
# 
#   }
# 
#   for (g in 2:heightTop)
#   {
# 
#     full_imgTop_dat = cbind(full_imgTop_dat, as.numeric(imgTop_dat[,,g]))
#   }
# 
#   norm_full_imgFront_dat = full_imgFront_dat/max(full_imgFront_dat)
#   norm_full_imgSide_dat = full_imgSide_dat/max(full_imgSide_dat)
#   norm_full_imgTop_dat = full_imgTop_dat/max(full_imgTop_dat)
# 
#   id = sprintf("%04d", j)
# 
#   listOfPixels[[paste0("image_",id,"_front")]] = norm_full_imgFront_dat
#   listOfPixels[[paste0("image_",id,"_side")]] = norm_full_imgSide_dat
#   listOfPixels[[paste0("image_",id,"_top")]] = norm_full_imgTop_dat
# }
# 
# saveRDS(listOfPixels, "pixelsMatrix")

data = readRDS("pixelsMatrix")

# reconstructed = Image(norm_full_img_dat, colormode = Grayscale)
# 
# display(reconstructed)
# 
# kmeans = kmeans(norm_full_img_dat, 10)
# 
# M = 200
# 
# height_seq = seq(0, height, length = M)
# width_seq = seq(0, width, length = M)
# height_rep = rep(height_seq, M)
# width_rep = rep(width_seq, each = M)
# 
# plot(height_rep ~ width_rep)
