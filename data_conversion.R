library(Cardinal)

# set wd to src file location
# get mass range & est. mass res.
# files <- list.files(path="./Normal", pattern="\\.imzML$", full.names=TRUE, recursive=FALSE)
# for (i in files) {
#   capture.output(f <- readMSIData(i),file='./out.txt',
#                  type = c("message"),
#                  append=TRUE)
# }



# load .imzML data and perform TIC normalization
read_imzML <- function(f) {
  f <- readMSIData(f, mass.range=c(160,1000),
                   resolution=100,
                   units="ppm")
  f <- normalize(f, method="tic")
  f <- smoothSignal(f, method="gaussian")
  f <- reduceBaseline(f, method="locmin")
  f <- process(f)
  f <- mzBin(f, from=160, to=1000, resolution=5000, units="ppm") %>% process()
  return(f)
}
  # f_peaks <- peakPick(f, method="mad")
  # f_peaks <- peakAlign(f_peaks, tolerance=200, units="ppm")
  # f_peaks <- peakFilter(f_peaks, freq.min=0.1)
  # f_peaks <- process(f_peaks)
  # f_peaks <- peakBin(f, ref=mz(f_peaks), type="height")
  # f_peaks <- f_peaks %>% process()
  # plot(f, coord=list(x=10, y=10))

# extract mass spectra and position information
extract <- function(f) {
  # f = dataframe from read_imzML
  # get coords for each pixel
  # XCoord <- coord(f)$x
  # YCoord <- coord(f)$y
  # X = max(XCoord)
  # Y = max(YCoord)
  # ZCoord <- coord(f)$z
  # coord <- data.frame(x=XCoord,y=YCoord,z=ZCoord)
  
  # mass spectra for each pixel
  mz <- as.data.frame(f@featureData)
  df <- as.data.frame(spectra(f))
  rownames(df) <- mz$mz
  
  #rm(f,coord,mz,XCoord,YCoord,ZCoord)
  rm(f)
  return(df)
}

# loop through all files and export as .csv
f_N <- list.files(path="./Normal", pattern="\\.imzML$", full.names=TRUE, recursive=FALSE)
f_T <- list.files(path="./Tumor", pattern="\\.imzML$", full.names=TRUE, recursive=FALSE)

export <- function(coord, df, name){
  # write.csv(coord, paste0(name, "coord", ".csv"), row.names=TRUE)
  write.csv(df, paste0(name,".csv"), row.names=TRUE)
}

# Normal
for(file in f_N) {
  f <- read_imzML(file)
  export(df, sub('\\.imzML$', '', file))
}

# Tumor
for(file in f_T) {
  f <- read_imzML(file)
  df <- extract(f)
  export(df, sub('\\.imzML$', '', file))
}

# extract (X, Y)
name_list = c(f_N, f_T)
X = c()
Y = c()
for(file in name_list) {
  f <- read_imzML(file)
  XCoord <- coord(f)$x
  YCoord <- coord(f)$y
  x = max(XCoord)
  y = max(YCoord)
  X <- append(X,x)
  Y <- append(Y,y)
}

coords <- data.frame(file = name_list, x = X, y = Y)
write.csv(coords, file = "coords.csv")

# t <- read_imzML("./Normal/202012122000-75_N_POS.imzML")
# dt <- extract(t)
# 
# XCoord <- coord(t)$x
# YCoord <- coord(t)$y
# ZCoord <- coord(t)$z
# coord <- data.frame(x=XCoord,y=YCoord,z=ZCoord)

# spectra(t)
