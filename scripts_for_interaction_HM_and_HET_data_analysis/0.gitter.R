#################################################
#           Read images of plates with gitter   #
#################################################

gitter.batch(image.files = "2018_07_S2", plate.format = c(32, 48),fast = c(4000), grid.save = "2018_07_S2/data", verbose = "p")
gitter.batch(image.files = "2018_08_MTX1_AM/J4", plate.format = c(32, 48),fast = c(4000), grid.save = "2018_08_MTX1_AM/data", verbose = "p")
gitter.batch(image.files = "2018_08_DMSO2_AM/J4", plate.format = c(32, 48),fast = c(4000), grid.save = "2018_08_MTX2_AM/data", verbose = "p")
