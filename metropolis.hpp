/*
 *
 *
 */

#include "vmmlib/vmmlib.hpp"
#include <stdlib.h>
#include <time.h>
#include <cmath>

#define DEBUG
#define ESPINES



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
  double probabilidad(double E_inicial, double E_final);
  double energia();
  double entropia();
  void print(void);
  int **obtener_primeros_vecinos(int fila, int columna);
  void set_temp(double new_value);
  void set_kb(double new_value);
  void set_A_prob(double new_value);
  void set_influencia_externa(double new_value);
  void set_condicion_externa(double new_value);
  void set_influencia_primeros_vecinos(double new_value);

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
  T valor_anterior = estado(fila,columna);
  double E_inicial = this->energia();
  this->estado(fila,columna) = this->posibles_estados[indice_estado];
  double E_final = this->energia();
  double prob = this->probabilidad(E_inicial, E_final);
  if ((rand()*1.0/RAND_MAX) >= prob) this->estado(fila,columna) = valor_anterior;
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
      T estado_ = this->estado(fila,columna);
      int **vecinos;      
      vecinos = obtener_primeros_vecinos(fila, columna);
#ifdef ESPINES
      E += this->influencia_externa*this->condicion_externa*estado_;
      for (int i=0; i<4; i++) 
      {
        T vecino_ = this->estado(vecinos[i][0], vecinos[i][1]);
        E += this->influencia_primeros_vecinos*vecino_/2; //Divido por 2 porque se están repitiendo emparejamientos
      }
#else 
      if ( == this->condicion_externa)
      {
        E -= this->influencia_externa;
      }
      for (int i=0; i<4; i++) 
      {
        T vecino_ = this->estado(vecinos[i][0], vecinos[i][1]);
        if (vecino_ == estado_) E -= this->influencia_primeros_vecinos;
      }
#endif //ESPINES
    }
  }
  return E;
}

TEMPLATE
double MODELO::probabilidad(double E_inicial, double E_final)
{
  double prob=-1.0;
  if (E_inicial>E_final) prob=this->A_prob;
  else{
    if (this->temp != 0 && this->kb != 0) prob=this->A_prob*exp(abs(E_inicial-E_final)/this->temp*this->kb);
    else prob=0;
  }
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
  return NULL;
}

TEMPLATE
int** MODELO::obtener_primeros_vecinos(int fila, int columna)
{
  int** vecinos = new int*[4];
  for (int vecino=0; vecino<4; vecino++)
  {
    vecinos[vecino] = new int[2];
  }
  //Vecino de arriba.
  vecinos[0][0] = (fila>0)? (fila-1):(this->filas-1);
  vecinos[0][1] = columna;
  //Vecino de la derecha.
  vecinos[1][0] = fila;
  vecinos[1][1] = (columna<(this->columnas-1))? (columna+1):0;
  //Vecino de abajo.
  vecinos[2][0] = (fila<(this->filas-1))? (fila+1):0;
  vecinos[2][1] = columna;
  //Vecino de la izquierda.
  vecinos[3][0] = fila;
  vecinos[3][1] = (columna>0)? (columna-1):(this->columnas-1);

  return vecinos;
}



//getters y setters
TEMPLATE
void MODELO::set_temp(double new_temp)
{
  this->temp = new_temp;
}

TEMPLATE
void MODELO::set_kb(double new_kb)
{
  this->kb = new_kb;
}

TEMPLATE
void MODELO::set_influencia_externa(double new_influencia_externa)
{
  this->influencia_externa = new_influencia_externa;
}

TEMPLATE
void MODELO::set_A_prob(double new_A)
{
  this->A_prob = new_A;
}

  double condicion_externa;
  double influencia_primeros_vecinos;

TEMPLATE
void MODELO::set_condicion_externa(double new_condicion_externa)
{
  this->condicion_externa = new_condicion_externa;
}

TEMPLATE
void MODELO::set_influencia_primeros_vecinos(double new_influencia_primeros_vecinos)
{
  this->influencia_primeros_vecinos = new_influencia_primeros_vecinos;
}

