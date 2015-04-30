/*
 *
 *
 *
 */

#include <iostream>
#include <cstdlib>
#include "metropolis.hpp"
#define DEF_NOMBRE(_A_) "graphs"+(char)_A_+".plot"


const int size_x = 10;
const int size_y = 10;
modelo<size_y,size_x,string> modelo1;
void set_ceros(double **matriz, int filas, int columnas);
void antisimetriza(double **matriz, int filas, int columnas);
void simulacion(int N_posibles_estados, string *posibles_estados, int numero_pasos, int longitud_paso);
void simulacion2(int N_posibles_estados, string *posibles_estados, int numero_pasos, int N_temp, double delta_temp);

int main(int argc, char** argv){
  
  double temp,
	A_prob,
	influencia_externa;
  string condicion_externa;

  const int N_temp = 1000;
  const   int delta_temp = 0;

  //parametros
  const int N_posibles_estados = 2;
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
/*
  posibles_estados[0] = "egoista";
  posibles_estados[1] = "altruista";
  posibles_estados[2] = "precavido";
*/
  posibles_estados[0] = "up";
  posibles_estados[1] = "down";

  modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
  modelo1.llenar();  
  if (argc==1){
	temp=0;
	A_prob=1;
	influencia_externa=0;
	condicion_externa=posibles_estados[0];	
  }else if (argc==2){
	temp=(double)atof(argv[1]);
	A_prob=1;
	influencia_externa=0;
	condicion_externa=posibles_estados[0];	
  }else{
	temp=(double)atof(argv[1]);
	A_prob=(double)atof(argv[2]);
	influencia_externa=(double)(atof(argv[3]));
	condicion_externa=argv[4];	
  }
  modelo1.set_temp(temp);
  modelo1.set_kb(1);
  modelo1.set_A_prob(A_prob);
  modelo1.set_influencia_externa(influencia_externa);
  modelo1.set_condicion_externa(condicion_externa);
/*
  influencia_vecinos[0][0] = -10;
  influencia_vecinos[1][1] = 5;
  influencia_vecinos[0][1] = 20;
  influencia_vecinos[0][2] = 0;
  influencia_vecinos[1][2] = 20;
  influencia_vecinos[2][2] = 0;
*/
  influencia_vecinos[0][0] = 1;
  influencia_vecinos[1][1] = 1;
  influencia_vecinos[0][1] = -1;
  antisimetriza(influencia_vecinos,N_posibles_estados,N_posibles_estados);
  const int longitud_paso=10;
  const int numero_pasos=250;
  modelo1.set_influencia_primeros_vecinos(influencia_vecinos);

  //simulacion(N_posibles_estados, posibles_estados, numero_pasos, longitud_paso);
  simulacion2(N_posibles_estados,posibles_estados, numero_pasos, N_temp, delta_temp);

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
	if (filas!=columnas) {std::cerr<<"La matriz tiene que ser cuadrada (antisimetrizacion).\n";std::exit(EXIT_FAILURE);}
  for (int fila=0; fila<filas; fila++)
  {
    for (int columna=fila+1; columna<columnas; columna++)
    {
      matriz[fila][columnas-columna] = -matriz[fila][columna];
    }
  }
}

void simulacion(int N_posibles_estados, string *posibles_estados, int numero_pasos, int longitud_paso)
  {
	 ofstream energia_t,
           entropia_t,
	   proporciones_t,
	   graphs;

  energia_t.open("energia_t.txt");
  entropia_t.open("entropia_t.txt");
  proporciones_t.open("proporciones_t.txt");
  graphs.open("graphs.plot");

	graphs<<"set terminal png size 640,300\nset output \"energia-beta"<<modelo1.get_temp()*modelo1.get_kb()<<"-A"<<modelo1.get_A_prob()<<"-inf"<<modelo1.get_influencia_externa()<<"-"<<modelo1.get_condicion_externa()<<".png\"\nset grid\nplot \"energia_t.txt\" with lines\nset output \"proporciones-beta"<<modelo1.get_temp()*modelo1.get_kb()<<"-A"<<modelo1.get_A_prob()<<"-inf"<<modelo1.get_influencia_externa()<<"-"<<modelo1.get_condicion_externa()<<".png\"\nset yrange [0:100]\nset grid\n plot ";

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
  for (int i=0; i<numero_pasos; i++) 
  {
	k++;
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
  }

  energia_t.close();
  entropia_t.close();
  proporciones_t.close();
  graphs.close();

  }


void simulacion2(int N_posibles_estados, string *posibles_estados, int numero_pasos, int N_temp, double delta_temp)
  {
	 ofstream magnetizacion,
	          graphs;

  magnetizacion.open("magnetizacion-T.txt");
  graphs.open("graphs.plot");
  
  magnetizacion<<"#Temperatura\tMagnetizacion";

	graphs<<"set terminal png size 640,300\nset output \"magnetizacion-beta"<<modelo1.get_temp()*modelo1.get_kb()<<"-A"<<modelo1.get_A_prob()<<"-inf"<<modelo1.get_influencia_externa()<<"-"<<modelo1.get_condicion_externa()<<".png\"\nset grid\nplot \"magnetizacion.txt\" with lines";


  cout<<"\n";
  int k=-1;
  int* num_estados;
  double temperatura;
  for (int temp_=0; temp_<N_temp; temp_++)
  {
	  for (int i=0; i<numero_pasos; i++) 
	  {
		  modelo1.cambiar_estado();
	  }
    num_estados = modelo1.get_cantidad_estado();
    temperatura = modelo1.get_temp();
    magnetizacion<<temperatura<<"\t"<<abs(num_estados[0]-num_estados[1]);
    modelo1.set_temp(temperatura+delta_temp*temp_);
  }
  magnetizacion.close();
  graphs.close();

  }
