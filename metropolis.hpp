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
  inline T& operator()( size_t row_index, size_t col_index );
  inline const T& operator()( size_t row_index, size_t col_index ) const;
  inline T& at( size_t row_index, size_t col_index );
  inline const T& at( size_t row_index, size_t col_index ) const;

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
  this->filas = (int)M;
  this->columnas = (int)N;
}


template< size_t M, size_t N, typename T>
void modelo<M,N,T>::definir_posibles_estados(int N_posibles_estados_in, double *posibles_estados_in)
{
  this->posibles_estados = posibles_estados_in;
  this->N_posibles_estados = N_posibles_estados_in;
}


template< size_t M, size_t N, typename T>
void modelo<M,N,T>::llenar()
{
  srand(time(NULL));
  srand(rand());
  int random_num=0;
  for (int fila=0; fila < this->filas; fila++)
  {
    for (int columna=0; columna < this->columnas; columna++)
    {
      random_num = rand()%(this->N_posibles_estados);
      this->estado(fila, columna) = this->posibles_estados[random_num];
    }
  }
}



template< size_t M, size_t N, typename T >
inline T&
modelo< M, N, T >::at( size_t row_index, size_t col_index )
{
    return this->estado(row_index, col_index);
}


template< size_t M, size_t N, typename T >
const inline T&
modelo< M, N, T >::at( size_t row_index, size_t col_index ) const
{
    return this->estado(row_index, col_index);
}

template< size_t M, size_t N, typename T >
inline T&
modelo< M, N, T >::operator()( size_t row_index, size_t col_index )
{
    return this->estado(row_index, col_index);
}


template< size_t M, size_t N, typename T >
const inline T&
modelo< M, N, T >::operator()( size_t row_index, size_t col_index ) const
{
    return this->estado(row_index, col_index);
}
