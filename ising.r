###
# Parámetros del programa
###


N_filas <- 20
N_columnas <- 20
J0 <- 10
H0 <- 0
T0 <- 0
kb <- 1e-10
set.seed(cat(Sys.time())) #cambia la semilla de los numeros aleatorios
posibles_estados <- c(-1/2, 1/2)
Sp <- matrix(nrow=N_filas, ncol=N_columnas)


###
# Funciones
###


adyacentes <- function(fila, columna){
  fila_adyacente <- 1
  columna_adyacente <- 1
  if (fila<N_filas){
    fila_adyacente <- fila+1
  }
  if (columna<N_columnas){
    columna_adyacente <- columna+1
  }
  return(c(fila_adyacente, columna_adyacente))
}

energia <- function(estado,termino_de_acoplo=J0,influencia_externa=H0) {
  energia_calc <- 0
  for (fila in 1:length(estado[,1])){
    for (columna in 1:length(estado[1,])){
      adyacentes_ <- adyacentes(fila, columna)
      fila_adyacente <- adyacentes_[1]
      columna_adyacente <- adyacentes_[2]
      energia_calc <- energia_calc - termino_de_acoplo*(estado[fila][columna]*estado[fila_adyacente][columna] + estado[fila][columna]*estado[fila][columna_adyacente]) - influencia_externa*estado[fila][columna]
    }
  }
  return(energia_calc)
}

cambiar_estado_aleatoriamente <- function(estado, posibles_estados_=posibles_estados){
  nuevo_estado <- estado
  fila_aleatoria <- sample(1:N_filas,1)
  columna_aleatoria <- sample(1:N_columnas,1)
  nuevo_estado[fila_aleatoria][columna_aleatoria] = sample(posibles_estados_,1,replace=TRUE)
  return(nuevo_estado)  
}

probabilidad <- function(estado, nuevo_estado, parametro_probabilidad=1, termino_de_acoplo=J0, influencia_externa=H0, influencia_del_ruido=kb, ruido=T0){
    energia_inicial <- energia(estado,termino_de_acoplo,influencia_externa)
    energia_final <- energia(nuevo_estado,termino_de_acoplo,influencia_externa)
    if (energia_final <= energia_inicial){
      return(parametro_probabilidad)
    }else{A 
      prob <- parametro_probabilidad*exp((energia_final-energia_inicial)/(influencia_del_ruido*ruido))
      return(prob)
    }
}

entropy <- function(A){
  return(0)
}

llenar <- function(Sp, posibles_estados, N_filas, N_columnas){
  for (fila in 1:N_filas){
    Sp[fila,] <<- sample(posibles_estados,N_columnas,replace=TRUE)
  }
}

tests <- function(){
  # Se rellena la matriz con estados aleatorios, escogidos entre los posibles.
  llenar(Sp,posibles_estados,N_filas,N_columnas)
  E <- energia(Sp)
  print(E)
  return(NULL)
}



###
# Parte principal del programa
###
main <- function(){
  # Se rellena la matriz con estados aleatorios, escogidos entre los posibles.
  llenar(Sp,posibles_estados,N_filas,N_columnas)
  
  estado <- Sp
  paso <- 0
  energias_libres <- matrix(nrow=1, ncol=N_filas*N_columnas)
  entropias <- matrix(nrow=1, ncol=N_filas*N_columnas)
  pasos <- matrix(nrow=1, ncol=N_filas*N_columnas)
  
  for (fila in 1:N_filas){
    for (columnas in 1:N_columnas){
      paso <- paso+1
      nuevo_estado <- cambiar_estado_aleatoriamente(Sp, posibles_estados)
      prob <- probabilidad(estado, nuevo_estado)
      if (runif(1) <= prob){
        estado <- nuevo_estado
      }
      energias_libres[paso] <- energia(estado)
      entropias[paso] <- entropia(estado)
      pasos[paso] <- paso
    }
  }
  
  plot(pasos, energias_libres)
}

tests()