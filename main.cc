/*
 *
 *
 *
 */

#include <iostream>
using namespace std;
#include "metropolis.hpp"

//#define ESTADOS2
//#define ESTADOS3
#define ESTADOS4
#define PROGRESIVO


void set_ceros(double **matriz, int filas, int columnas);
void antisimetriza(double **matriz, int filas, int columnas);

int main(int argc, char** argv){
  double temp=0;
  if (argc==2) temp=(double)atof(argv[1]);
  //parametros
#ifdef ESTADOS2
  const int N_posibles_estados = 2;
#endif
#ifdef ESTADOS3
  const int N_posibles_estados = 3;
#endif
#ifdef ESTADOS4
  const int N_posibles_estados = 4;
#endif
  string posibles_estados[N_posibles_estados];
  double **influencia_vecinos = new double*[N_posibles_estados];
  for (int i=0; i<N_posibles_estados; i++)
  {
    influencia_vecinos[i] = new double[N_posibles_estados];
    for (int j=0; j<N_posibles_estados; j++)
    {
      influencia_vecinos[i][j] = 0;
    }
  }
  const int size_x = 100;
  const int size_y = 100;
  modelo<size_y,size_x,string> modelo1;
#ifdef ESTADOS2
  posibles_estados[0] = "altruista";
  posibles_estados[1] = "egoista";
#endif
#ifdef ESTADOS3
  posibles_estados[0] = "altruista";
  posibles_estados[1] = "egoista";
  posibles_estados[2] = "copia";
#endif
#ifdef ESTADOS4
  posibles_estados[0] = "altruista";
  posibles_estados[1] = "egoista";
  posibles_estados[2] = "copia";
  posibles_estados[3] = "contrario";
#endif
  modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
  modelo1.llenar();
  modelo1.set_temp(temp);
  modelo1.set_kb(1);
  modelo1.set_A_prob(1);
  modelo1.set_influencia_externa(5);
#ifdef ESTADOS2
  modelo1.set_condicion_externa("altruista");
  influencia_vecinos[0][0] = 10;
  influencia_vecinos[1][1] = -5;
  influencia_vecinos[0][1] = -20;
#endif
#ifdef ESTADOS3
  modelo1.set_condicion_externa("altruista");
  influencia_vecinos[0][0] = 10;
  influencia_vecinos[1][1] = -5;
  influencia_vecinos[2][2] = 0;
  influencia_vecinos[0][1] = -20;
  influencia_vecinos[0][2] = 10;
  influencia_vecinos[1][2] = -5;
#endif
#ifdef ESTADOS4
  modelo1.set_condicion_externa("copia");
  influencia_vecinos[0][0] = 10;
  influencia_vecinos[1][1] = -5;
  influencia_vecinos[2][2] = 0;
  influencia_vecinos[3][3] = 0;
  influencia_vecinos[0][1] = -20;
  influencia_vecinos[0][2] = 10;
  influencia_vecinos[1][2] = -5;
  influencia_vecinos[0][3] = -20;
  influencia_vecinos[1][3] = 20;
  influencia_vecinos[2][3] = 0;
#endif
  antisimetriza(influencia_vecinos,N_posibles_estados,N_posibles_estados);
  const int longitud_paso=10;
  const int numero_pasos=1000;
  modelo1.set_influencia_primeros_vecinos(influencia_vecinos);




  ofstream energia_t,
           entropia_t,
	   proporciones_t,
	   graphs;
  energia_t.open("energia_t.txt");
  entropia_t.open("entropia_t.txt");
  proporciones_t.open("proporciones_t.txt");
  graphs.open("graphs.plot");

	graphs<<"set terminal png size 640,300\nset output \"energia-kbT"<<modelo1.get_temp()*modelo1.get_kb()<<"-A"<<modelo1.get_A_prob()<<"-inf"<<modelo1.get_influencia_externa()<<"-"<<modelo1.get_condicion_externa()<<".png\"\nset grid\nplot \"energia_t.txt\" with lines\nset output \"proporciones-kbT"<<modelo1.get_temp()*modelo1.get_kb()<<"-A"<<modelo1.get_A_prob()<<"-inf"<<modelo1.get_influencia_externa()<<"-"<<modelo1.get_condicion_externa()<<".png\"\nset yrange [0:100]\nset grid\n plot ";

  energia_t<<"#paso\tenergia\n";
  entropia_t<<"#paso\tentropia\n";
  proporciones_t<<"#paso\t";
  for (int estado=0; estado<N_posibles_estados; estado++)
  {
    graphs<<"\"proporciones_t.txt\" using 1:"<<estado+2<<" title \""<<posibles_estados[estado]<<"\" with lines";
    if (estado<N_posibles_estados-1) graphs<<", ";
    proporciones_t<<posibles_estados[estado]<<"\t";
  }
  proporciones_t<<"\n";

  int* num_estados;
  cout<<"\n";
  int k=-1;
#ifdef PROGRESIVO
  int hk = -1;
#endif
  for (int i=0; i<numero_pasos; i++) 
  {
	k++;
#ifdef PROGRESIVO
        hk++;
#endif
	modelo1.cambiar_estado();
        if (k==longitud_paso){
		modelo1.contar_estados();
		num_estados = modelo1.get_cantidad_estado();
          	proporciones_t<<i<<"\t";
		for (int j=0; j<N_posibles_estados; j++)
		{
			proporciones_t<<num_estados[j]*100.0/(size_x*size_y)<<"\t";
		}
		proporciones_t<<"\n";
		energia_t<<i<<"\t"<<modelo1.energia()<<"\n";
                k=0;
	}  
	if (hk==100){
#ifdef PROGRESIVO
	  influencia_vecinos[0][0] = rand()%50-100;
	  influencia_vecinos[1][1] = rand()%50-100;
	  influencia_vecinos[2][2] = rand()%50-100;
	  influencia_vecinos[3][3] = rand()%50-100;
	  influencia_vecinos[0][1] = rand()%50-100;
	  influencia_vecinos[0][2] = rand()%50-100;
	  influencia_vecinos[1][2] = rand()%50-100;
	  influencia_vecinos[0][3] = rand()%50-100;
	  influencia_vecinos[1][3] = rand()%50-100;
	  influencia_vecinos[2][3] = rand()%50-100;
	  antisimetriza(influencia_vecinos,N_posibles_estados,N_posibles_estados);
	  modelo1.set_influencia_primeros_vecinos(influencia_vecinos);
#endif
	}
  }

  energia_t.close();
  entropia_t.close();
  proporciones_t.close();
  graphs.close();
}


void set_ceros(double **matriz, int filas, int columnas)
{
  for (int fila=0; fila<filas; fila++)
  {
    for (int columna=0; columna<columnas; columna++)
    {
      matriz[fila][columna] = 0;
    }
  }
}

void antisimetriza(double **matriz, int filas, int columnas)
{
  if (filas!=columnas) {std::cerr<<"La matriz tiene que ser cuadrada (antisimetrizacion).\n";}
  for (int fila=0; fila<filas; fila++)
  {
    for (int columna=fila+1; columna<columnas; columna++)
    {
      matriz[fila][columnas-columna] = -matriz[fila][columna];
    }
  }
}
