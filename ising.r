# runif (Number of elements, start, end) FLOAT
# sample (start:end, number of elements) INT or ELEMENT IN ARRAY

Nx <- 20
Ny <- 20
J0 <- 10
H0 <- 0


energy <- function(A,J=J0,H=H0) {
  E <- 0
  for (i in 1:length(A[,1])){
    for (j in 1:length(A[1,])){
      ii <- 1
      jj <- 1
      if (i<Ny){
        ii <- i+1
      }
      if (j<Nx){
        jj <- j+1
      }
      E <- E - J*(A[i,j]*A[ii,j] + A[i,j]*A[i,jj]) - H*A[i,j]
    }
  }
  return(E)
}

entropy <- function(A){
  
}




Sp <- matrix(nrow=Nx, ncol=Ny)
# Se rellena la matriz con +0.5 o -0.5 aleatoriamente
for (i in 1:Ny){
  for (j in 1:Nx){
    Sp[i,j] <- sample(c(-0.5,0.5),1)
  }
}



t <- rep(NA, Nx*Ny)
En <- rep(NA, Nx*Ny)

for (i in 1:length(Sp[,1])){
  for (j in 1:length(Sp[1,])){
    E1 <- energy(Sp)
    Sp[i,j] <- Sp[i,j] * -1
    E2 <- energy(Sp)
    E <- E2
    if (E1 < E2){
      Sp[i,j] <- Sp[i,j] * -1
      E <- E1
    }
    En[j+Ny*(i-1)] <- E
    t[j+Ny*(i-1)] <- j+Ny*(i-1)
  }
}

plot(t,En,"l")