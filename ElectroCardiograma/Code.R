rm(list=ls())

###########################
### Leer el archivo hea ###
###########################

hea_file <- "Data/ECGPCG0003.hea"
hea_content <- readLines(hea_file)

# Mostrar contenido del archivo .hea
print(hea_content)


#Extrear datos del archivo .hea
n_muestra<-240000 # Muestra
freq<-8000 # Frecuenca (cada momento que se toma una medida)
gain_ecg <- 110554.8863 # las ganancias

num_signals <- 2  #el número de señales en el archivo (ECG y PCG)


###########################
### Leer el archivo dat ###
###########################



dat_file <- "Data/ECGPCG0003.dat"


# Leer los datos como enteros de 16 bits (big-endian si es necesario)
raw_data <- readBin(dat_file, 
                    what = "integer", 
                    size = 2, 
                    n=n_muestra*num_signals,
                    signed = TRUE, 
                    endian = "little")


# Reestructurar los datos en dos columnas (una por cada señal: ECG y PCG)
data_matrix <- matrix(raw_data, ncol = num_signals, byrow = TRUE)

ECG<-data_matrix[,1]
length(ECG)

##########################
## Escala de conversión ##
##########################



ECG_mV <- ECG/gain_ecg

ECG_data <- data.frame(Time = seq(0, length(ECG_mV) - 1) / freq,  # Tiempo en segundos
                       ECG = ECG_mV)

#################### 
##### Graficar ####
####################

library(ggplot2)

# Graficar la señal ECG
ggplot(ECG_data, aes(x = Time, y = ECG)) +
  geom_line(color = "blue") +
  labs(title = "Señal de ECG", x = "Tiempo (s)", y = "Amplitud (mV)") +
  theme_minimal()


#################################
####### DETECTAR PICOS R ########
#################################

install.packages("pracma")  
library(pracma)
# Definir un umbral basado en la señal ECG (puedes ajustarlo según la amplitud de tu señal)
umbral <- max(ECG_mV) * 0.5  # La mitad del pico más alto

# Detectar picos R en la señal ECG
picos <- findpeaks(ECG_mV, 
                   minpeakheight = umbral, 
                   minpeakdistance = freq * 0.6)


# Extraer tiempos de los picos
tiempo_picos <- ECG_data$Time[picos[, 2]]
amplitud_picos <- picos[, 1]

# Ordenarlos 
orden_picos <- order(tiempo_picos)
tiempo_picos_ordenados <- tiempo_picos[orden_picos]
amplitud_picos_ordenados <- amplitud_picos[orden_picos]


# Crear un data frame con los picos detectados
picos_df <- data.frame(Time = tiempo_picos_ordenados, 
                       Amplitud = amplitud_picos_ordenados)

# Graficar la señal ECG con los picos R resaltados
ggplot(ECG_data, aes(x = Time, y = ECG)) +
  geom_line(color = "blue") +  # Señal ECG
  geom_point(data = picos_df, aes(x = Time, y = Amplitud), color = "red", size = 2) +  # Picos R
  labs(title = "Señal ECG con detección de picos R", x = "Tiempo (s)", y = "Amplitud (mV)") +
  theme_minimal()

##########################################
####### CALCULA LOS INTERVALOS RR ########
##########################################



# Calcular los intervalos RR (diferencias entre los tiempos de los picos R)
intervalos_RR <- diff(tiempo_picos_ordenados)

# Estadísticas básicas de los intervalos RR
rr_mean <- mean(intervalos_RR)  # Promedio
rr_sd <- sd(intervalos_RR)      # Desviación estándar
rr_min <- min(intervalos_RR)    # Mínimo
rr_max <- max(intervalos_RR)    # Máximo

# Mostrar resultados
cat("Intervalo RR promedio:", rr_mean, "s\n")
cat("Desviación estándar:", rr_sd, "s\n")
cat("Intervalo RR mínimo:", rr_min, "s\n")
cat("Intervalo RR máximo:", rr_max, "s\n")

# Crear un dataframe para graficar
RR_data <- data.frame(Time = tiempo_picos[-1], RR_Interval = intervalos_RR)

# Graficar los intervalos RR
ggplot(RR_data, aes(x = Time, y = RR_Interval)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "red", size = 1) +
  labs(title = "Intervalos RR a lo largo del tiempo", x = "Tiempo (s)", y = "Intervalo RR (s)") +
  theme_minimal()

#####################################
##### DETECTAR VALORES ANOMALOS #####
#####################################


# Calcular los cuartiles
Q1 <- quantile(ECG_mV, 0.25)
Q3 <- quantile(ECG_mV, 0.75)
IQR <- Q3 - Q1

# Definir los límites inferior y superior para detectar outliers
limite_inferior <- Q1 - 1.5 * IQR
limite_superior <- Q3 + 1.5 * IQR

# Detectar los valores anómalos
outliers_IQR <- which(ECG_mV < limite_inferior | ECG_mV > limite_superior)

# Mostrar los valores anómalos detectados
outliers_IQR_values <- ECG_mV[outliers_IQR]
outliers_IQR_times <- ECG_data$Time[outliers_IQR]

# Mostrar los valores y tiempos de los outliers
data.frame(Time = outliers_IQR_times, ECG_mV = outliers_IQR_values)

# Graficar la señal ECG con los outliers detectados usando IQR
ggplot(ECG_data, aes(x = Time, y = ECG_mV)) +
  geom_line(color = "blue") +  # Señal ECG
  geom_point(data = data.frame(Time = outliers_IQR_times, ECG_mV = outliers_IQR_values),
             aes(x = Time, y = ECG_mV), color = "red", size = 2) +  # Outliers en rojo
  labs(title = "Señal ECG con Outliers detectados (IQR)", x = "Tiempo (s)", y = "Amplitud (mV)") +
  theme_minimal()






