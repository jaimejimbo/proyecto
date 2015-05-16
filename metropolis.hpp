/*
 *  - Clase principal de la simulación. En ella se definen los estados y métodos 
 * principales.
 *
 */

#include <cstdlib>
#include <time.h>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
//#include <assert>
using namespace std;
//Depuración
#define NDEBUG
#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif
//#define ESPINES

#define LONG


//uso un template(plantilla) para que se puedan usar estados definidos como texto, números u objetos abstractos.
template< size_t M, size_t N, typename T = string>
class modelo
{
/*
 *	Clase modelo
 *	tiene los métodos y variables necesarias para realizar las simulaciones
 *     de una forma más cómoda.
*/


public:
  modelo();
  void definir_posibles_estados(int N_posibles_estados, T *posibles_estados);	
  //llena la matriz con estados escogidos al azar entre los posibles.
  void llenar();

  //Métodos de acceso a los valores de la matriz (se accede como matriz(fila,columna) en vez de matriz[fila][columna] empezando por (0,0).
  inline T& operator()( size_t row_index, size_t col_index );
  inline const T& operator()( size_t row_index, size_t col_index ) const;
  inline T& at( size_t row_index, size_t col_index );
  inline const T& at( size_t row_index, size_t col_index ) const;


  //métodos importantes
  //-------------------
  void cambiar_estado();
  //calcula la probabilidad con las fórmulas del algoritmo de Metrópolis.
  double probabilidad(double E_inicial, double E_final);
  double energia();
  //este método no está implementado
  double entropia();
  //convierte la matriz a una cadena de caracteres para mostrarla en pantalla con cout.
  string to_string();
  //devuelve un vector con los primeros vecinos (arriba, derecha, abajo, izquierda)
  int **obtener_primeros_vecinos(int fila, int columna);
  //mira cuantos elementos hay en cada estado y lo guarda en la variable cantidades_por_estado
  void contar_estados();


  //métodos para cambiar los valores de los parámetros, y para obtenerlos.
  //-------------------
  void set_temp(double new_value);
  void set_kb(double new_value);
  void set_A_prob(double new_value);
  void set_influencia_externa(double new_value);
  void set_condicion_externa(T new_value);
  void set_influencia_primeros_vecinos(double **new_value);
  double get_temp();
  double get_kb();
  double get_A_prob();
  double get_influencia_externa();
  T get_condicion_externa();
  int* get_cantidades_por_estado();


private:
  //parámetros del sistema (Temperatura, constante de Boltzman, parámetro A de las ecuaciones de Metrópolis.
  double temp;
  double kb;
  double A_prob;
  T **estado;
  //dimensiones de la matriz 2D
  int filas;
  int columnas;
  int N_posibles_estados;
  //posibles estados en los que puede estar un elemento del sistema
  T *posibles_estados;
  //Campo magnético externo o similar (y si influencia sobre el sistema): H,m
  double influencia_externa;
  T condicion_externa;
  //Como interactúan los elementos entre sí. Es una matriz (elemento1 con elemento1, elemento1 con elemento2...)
  double **influencia_primeros_vecinos;
  //cuántos estados hay de cada tipo
  int *cantidades_por_estado;
};











#define TEMPLATE                  template< size_t M, size_t N, typename T>
#define MODELO                    modelo<M,N,T>
#define ADYACENTE(pos, length)    pos<length?pos+1:0

TEMPLATE
MODELO::modelo()
{
/*
*	Creador
*/
  this->filas = (int)M;
  this->columnas = (int)N;
  estado = new T*[this->filas];
  for (int fila=0; fila<this->filas; fila++)
  {
    estado[fila] = new T[this->columnas];
  }
}


TEMPLATE
void MODELO::definir_posibles_estados(int N_posibles_estados_in, T *posibles_estados_in)
{
/*
*   Se definen los posibles estados.
*/
  this->posibles_estados = posibles_estados_in;
  this->N_posibles_estados = N_posibles_estados_in;
  this->cantidades_por_estado = new int[N_posibles_estados_in];
}

TEMPLATE
void MODELO::llenar()
{
/*
*  Llena la matriz de estados con estados aleatorios escogidos entre los posibles estados.
*/
  srand(time(NULL));
  srand(rand());
  int random_num=0;
  for (int fila=0; fila < this->filas; fila++)
  {
    for (int columna=0; columna < this->columnas; columna++)
    {
      random_num = rand()%(this->N_posibles_estados);
      this->estado[fila][columna] = this->posibles_estados[random_num];
    }
  }
}



TEMPLATE
inline T&
MODELO::at( size_t row_index, size_t col_index )
{
/*
*    Este método no lo uso para nada, pero he visto que en librerías de matrices lo definen para usarlo en el operador()
*/
    ASSERT(0<=row_index,"");
    ASSERT(row_index<this->filas,"");
    ASSERT(0<=col_index,"");
    ASSERT(col_index<this->columnas,"");
    return this->estado(row_index, col_index);
}


TEMPLATE
const inline T&
MODELO::at( size_t row_index, size_t col_index ) const
{
/*
*	Lo mismo que arriba
*/
    ASSERT(0<=row_index,"");
    ASSERT(row_index<this->filas,"");
    ASSERT(0<=col_index,"");
    ASSERT(col_index<this->columnas,"");
    return this->estado(row_index, col_index);
}

TEMPLATE
inline T&
MODELO::operator()( size_t row_index, size_t col_index )
{
/*
*   Permite acceder a los valores de la matriz con paréntesis.
*/
    ASSERT(0<=row_index,"");
    ASSERT(row_index<this->filas,"");
    ASSERT(0<=col_index,"");
    ASSERT(col_index<this->columnas,"");
    return this->estado(row_index, col_index);
}


TEMPLATE
const inline T&
MODELO::operator()( size_t row_index, size_t col_index ) const
{
/*
*	Lo mismo. A esta función se la llama cuando se usa en una zona de constantes.
*/
    ASSERT(0<=row_index,"");
    ASSERT(row_index<this->filas,"");
    ASSERT(0<=col_index,"");
    ASSERT(col_index<this->columnas,"");
    return this->estado(row_index, col_index);
}


TEMPLATE
void MODELO::cambiar_estado()
{
/*
*     Cambia el estado de una posición aleatoria por un estado aleatorio. Para ello usa la probabilidad(por lo que tiene que obtener las energías).
*/
  srand(time(NULL));
  srand(rand());
#ifdef LONG
  for (int ii=0; ii<this->filas*this->columnas; ii++)
  {
#endif
    int fila;
    int columna;
    double E_inicial = 0,
           E_final=0;
    int indice,
        indice_v;
	  E_inicial=0;
	  E_final=0;
    fila = (int)(rand()%(this->filas));
    columna = (int)(rand()%(this->columnas));
	  int indice_estado = (int)(rand()%(this->N_posibles_estados));
	  T valor_anterior = estado[fila][columna];
	  T nuevo_estado = this->posibles_estados[indice_estado];


/*	  indice = -1;
	  indice_v = -1;
	  T estado_ = this->estado[fila][columna];
	  for (int i=0; i<this->N_posibles_estados; i++){
	    if (this->posibles_estados[i] == estado_) indice=i;
	  }
	  int **vecinos;      
	  vecinos = obtener_primeros_vecinos(fila, columna);
	  if (estado_ == this->condicion_externa)
	  {
	    E_inicial -= this->influencia_externa;
	  }
	  for (int i=0; i<4; i++) 
	  {
	    T vecino_ = this->estado[vecinos[i][0]][vecinos[i][1]];
	    for (int j=0; j<this->N_posibles_estados; j++){
	      if (this->posibles_estados[j] == vecino_) indice_v=j;
	    }
	    ASSERT(indice>=0,"indice < 0");
	    ASSERT(indice<this->N_posibles_estados,"indice > N_posibles_estados");
	    ASSERT(indice_v>=0,"indice_v < 0");
	    ASSERT(indice_v<this->N_posibles_estados,"indice_v > N_posibles_estados");

	    E_inicial -= this->influencia_primeros_vecinos[indice][indice_v];
	  }
*/
    E_inicial = this->energia();
	  this->estado[fila][columna] = nuevo_estado;
    E_final = this->energia();
/*
	  indice = -1;
	  indice_v = -1;
	  estado_ = nuevo_estado;
	  for (int i=0; i<this->N_posibles_estados; i++){
	    if (this->posibles_estados[i] == estado_) indice=i;
	  }
	  if (nuevo_estado == this->condicion_externa)
	  {
	    E_final -= this->influencia_externa;
	  }
	  for (int i=0; i<4; i++) 
	  {
	    T vecino_ = this->estado[vecinos[i][0]][vecinos[i][1]];
	    for (int j=0; j<this->N_posibles_estados; j++){
	      if (this->posibles_estados[j] == vecino_) indice_v=j;
	    }
	    ASSERT(indice>=0,"indice < 0");
	    ASSERT(indice<this->N_posibles_estados,"indice > N_posibles_estados");
	    ASSERT(indice_v>=0,"indice_v < 0");
	    ASSERT(indice_v<this->N_posibles_estados,"indice_v > N_posibles_estados");

	    E_final -= this->influencia_primeros_vecinos[indice][indice_v];
	  }
*/  
	  //cout<<E_inicial<<"\t"<<E_final<<"\n";
	  double prob = this->probabilidad(E_inicial, E_final);
	  if ((rand()*1.0/RAND_MAX) >= prob) {this->estado[fila][columna] = valor_anterior;}
#ifdef LONG
  }
#endif
}

TEMPLATE
string MODELO::to_string()
{
/*
*    Permite mostrar la matriz en pantalla (para depurar)
*/
  string output="";
  for (int fila=0; fila<this->filas; fila++)
  {
    for (int columna=0; columna<this->columnas; columna++)
    {
      output += (string)this->estado[fila][columna]+"\t";
    }
    output+="\n";
  }
  return output;
}

TEMPLATE
double MODELO::energia()
{
/*
*	Calcula la energía mediante una expresión tipo Landau.
*/
  int indice, indice_v, fila, columna;
  double E=0;
  for (int fila=0; fila<this->filas; fila++)
  {
    for (int columna=0; columna<this->columnas; columna++)
    {
      indice = -1;
      indice_v = -1;
      T estado_ = this->estado[fila][columna];
      for (int i=0; i<this->N_posibles_estados; i++){
        if (this->posibles_estados[i] == estado_) indice=i;
      }
      int **vecinos;      
      vecinos = obtener_primeros_vecinos(fila, columna);
      if (estado_ == this->condicion_externa)    //no lo hago con sumatorios, sino con condicionales para poder generalizar el proceso
      {
        E -= this->influencia_externa;
      }
      for (int i=0; i<4; i++) 
      {
        T vecino_ = this->estado[vecinos[i][0]][vecinos[i][1]];
	for (int j=0; j<this->N_posibles_estados; j++){
	  if (this->posibles_estados[j] == vecino_) indice_v=j;
	}
	ASSERT(indice>=0,"indice < 0");
	ASSERT(indice<this->N_posibles_estados,"indice > N_posibles_estados");
	ASSERT(indice_v>=0,"indice_v < 0");
	ASSERT(indice_v<this->N_posibles_estados,"indice_v > N_posibles_estados");

        E -= this->influencia_primeros_vecinos[indice][indice_v];
      }
    }
  }
  return E;
}

TEMPLATE
double MODELO::probabilidad(double E_inicial, double E_final)
{
/*
*     Se calcula la probabilidad con el método de Metrópolis
*/
  double prob;
  if (E_inicial>=E_final) prob=this->A_prob;
  else{
    if (this->temp != 0 && this->kb != 0) prob=this->A_prob*exp(-abs(E_inicial-E_final)/(this->temp*this->kb));
    else prob=0.0;
  }
  ASSERT(0<=prob,"Fallo al calcular la probabilidad <0");
  ASSERT(prob<=1,"Fallo al calcular la probabilidad >1");
  return prob;
}

TEMPLATE
double MODELO::entropia()
{
/*
*	No implementado
*/
  return 0;
}

TEMPLATE
int** MODELO::obtener_primeros_vecinos(int fila, int columna)
{
/*
*	devuelve un vector con los primeros vecinos (^, ->, v, <-)
*/
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

TEMPLATE
void MODELO::contar_estados()
{
/*
*	Mira cuantos elementos están en cada estado.
*/
  for (int num_estado=0; num_estado<this->N_posibles_estados; num_estado++)
  {
    cantidades_por_estado[num_estado] = 0;  
  }
  for (int fila=0; fila<this->filas; fila++)
  {
    for (int columna=0; columna<this->columnas; columna++)
    {
      for (int num_estado=0; num_estado<this->N_posibles_estados; num_estado++)
      {
        T estado = this->estado[fila][columna];
        if (estado == posibles_estados[num_estado]) cantidades_por_estado[num_estado]++;
      }
    
    }
  }
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
  ASSERT(0<=new_A,"El parámetro A tiene que estar entre 0 y 1");
  ASSERT(new_A<=1,"El parámetro A tiene que estar entre 0 y 1");
  this->A_prob = new_A;
}

TEMPLATE
void MODELO::set_condicion_externa(T new_condicion_externa)
{
  this->condicion_externa = new_condicion_externa;
}

TEMPLATE
void MODELO::set_influencia_primeros_vecinos(double **new_influencia_primeros_vecinos)
{
  this->influencia_primeros_vecinos = new_influencia_primeros_vecinos;
}


TEMPLATE
int* MODELO::get_cantidades_por_estado()
{
  return this->cantidades_por_estado;
}

TEMPLATE
double MODELO::get_temp()
{
	return this->temp;
}

TEMPLATE
double MODELO::get_kb()
{
	return this->kb;
}

TEMPLATE
double MODELO::get_A_prob()
{
	return this->A_prob;
}

TEMPLATE
double MODELO::get_influencia_externa()
{
	return this->influencia_externa;
}

TEMPLATE
T MODELO::get_condicion_externa()
{
	return this->condicion_externa;
}
