/*
 *
 *
 */

#include "vmmlib/vmmlib.hpp"
#include <stdlib.h>
#include <time.h>


template< size_t M, size_t N, typename T = float>
class modelo
{
public:
  modelo();
  void definir_posibles_estados(int N_posibles_estados, double *posibles_estados);
  void llenar();



private:
  vmml::matrix<M,N,T> estado;
  int filas;
  int columnas;
  int N_posibles_estados;
  double *posibles_estados;
};




template< size_t M, size_t N, typename T>
modelo<M,N,T>::modelo()
{
  this->filas = M;
  this->columnas = N;
  this->llenar();
}

template< size_t M, size_t N, typename T>
void modelo<M,N,T>::llenar()
{
  srand(time(NULL));
  for (int fila=0; fila < this->filas; fila++)
  {
    for (int columna=0; columna < this->columnas; columna++)
    {
      int random_num = rand()%N_posibles_estados;
      this->estado(fila,columna) = posibles_estados[random_num];
    }
  }
}

template< size_t M, size_t N, typename T>
void modelo<M,N,T>::definir_posibles_estados(int N_posibles_estados, double *posibles_estados)
{
  this->posibles_estados = posibles_estados;
  this->N_posibles_estados = N_posibles_estados;
}
