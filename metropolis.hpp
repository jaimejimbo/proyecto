/*
 *
 *
 */

#include "vmmlib/vmmlib.hpp"
#include <stdlib.h>
#include <time.h>
#include <cmath>

#define DEBUG



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
  void cambiar_estado();
  double probabilidad(double E_inicial, double E_final, double influencia_ruido, double ruido, double param_prob);
  double energia();
  double entropia();
  void print(void);
  int **obtener_vecinos(int fila, int columna);

private:
  double temp;
  double kb;
  double A_prob;
  vmml::matrix<M,N,T> estado;
  int filas;
  int columnas;
  int N_posibles_estados;
  double *posibles_estados;
  double influencia_externa;
  double condicion_externa;
  double influencia_primeros_vecinos;
};











#define TEMPLATE                  template< size_t M, size_t N, typename T>
#define MODELO                    modelo<M,N,T>
#define ADYACENTE(pos, length)    pos<length?pos+1:0

TEMPLATE
MODELO::modelo()
{
  this->filas = (int)M;
  this->columnas = (int)N;
}


TEMPLATE
void MODELO::definir_posibles_estados(int N_posibles_estados_in, double *posibles_estados_in)
{
  this->posibles_estados = posibles_estados_in;
  this->N_posibles_estados = N_posibles_estados_in;
}


TEMPLATE
void MODELO::llenar()
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



TEMPLATE
inline T&
MODELO::at( size_t row_index, size_t col_index )
{
    return this->estado(row_index, col_index);
}


TEMPLATE
const inline T&
MODELO::at( size_t row_index, size_t col_index ) const
{
    return this->estado(row_index, col_index);
}

TEMPLATE
inline T&
MODELO::operator()( size_t row_index, size_t col_index )
{
    return this->estado(row_index, col_index);
}


TEMPLATE
const inline T&
MODELO::operator()( size_t row_index, size_t col_index ) const
{
    return this->estado(row_index, col_index);
}


TEMPLATE
void MODELO::cambiar_estado()
{
  srand(time(NULL));
  srand(rand());
  int fila = (int)(rand()%(this->filas));
  int columna = (int)(rand()%(this->columnas));
  int indice_estado = (int)(rand()%(this->N_posibles_estados));
  this->estado(fila,columna) = this->posibles_estados[indice_estado];
}

TEMPLATE
void MODELO::print()
{
  for (int fila=0; fila<this->filas; fila++)
  {
    for (int columna=0; columna<this->columnas; columna++)
    {
      cout<<this->estado(fila,columna)<<"\t";
    }
    cout<<"\n";
  }
}

TEMPLATE
double MODELO::energia()
{
  double E=0;
  for (int fila=0; fila<this->filas; fila++)
  {
    for (int columna=0; columna<this->columnas; columna++)
    {
      int **vecinos = new int[4][2];
      vecinos = obtener_vecinos(fila, columna);
      for (int vecino=0; vecino<4; vecino++)
      {
        E += this->influencia_externa*this->condicion_externa + this->influencia_primeros_vecinos*this->estado(vecinos[vecino][0],vecinos[vecino][1]);
      }
    }
  }
}

TEMPLATE
double MODELO::probabilidad(double E_inicial, double E_final, double influencia_ruido, double ruido, double param_prob)
{
  double prob=-1.0;
  if (E_inicial>E_final) prob=param_prob;
  else prob=param_prob*exp(abs(E_inicial-E_final)/ruido*influencia_ruido);
#ifdef DEBUG
  if (prob>=0 and prob<=1) return prob;
  else cerr<<"Algo falla al calcular la probabilidad";
#else
  return prob;
#endif //DEBUG
}

TEMPLATE
double MODELO::entropia()
{

}
