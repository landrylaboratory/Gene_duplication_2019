#################################################
#           Read images of plates with gitter   #
#################################################

##Images of PCA plates were analysed using gitter 
#(Wagih and Parts, 2014) to quantify colony sizes. 
#(Contact authors for PCA pictures data)

rm(list=ls())

require(gitter)

#set the working directory to dir
setwd("")

#Run Gitter on image files
gitter.batch(image.files = "2018_07_S2", plate.format = c(32, 48),fast = c(4000), grid.save = "2018_07_S2/data", verbose = "p")
gitter.batch(image.files = "2018_08_MTX1_AM/J4", plate.format = c(32, 48),fast = c(4000), grid.save = "2018_08_MTX1_AM/data", verbose = "p")
gitter.batch(image.files = "2018_08_DMSO2_AM/J4", plate.format = c(32, 48),fast = c(4000), grid.save = "2018_08_MTX2_AM/data", verbose = "p")
