archivo_hea <- "ECGPCG0003.hea"
datos <- read.table(archivo_hea, 
                    header = FALSE, 
                    sep = "\t")
file_dat <- readBin("ECGPCG0003.dat", 
                      what = "integer", 
                      size = 2, 
                      n = 240000*2, 
                      endian = "little")
library(R.matlab)




head(datos)
