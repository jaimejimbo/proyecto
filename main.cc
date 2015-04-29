/*
 *
 *
 *
 */

#include <iostream>
using namespace std;
#include "metropolis.hpp"

void set_ceros(double **matriz, int filas, int columnas);

int main(int argc, char** argv){
  const int N_posibles_estados = 3;
  string posibles_estados[N_posibles_estados];
  posibles_estados[0] = "omnivoro";
  posibles_estados[1] = "carnivoro";
  posibles_estados[2] = "vegetariano";
  double **influencia_vecinos = new double*[3];
  for (int i=0; i<3; i++)
  {
    influencia_vecinos[i] = new double[3];
    for (int j=0; j<3; j++)
    {
      influencia_vecinos[i][j] = 0;
    }
  }
  set_ceros(influencia_vecinos,3,3);



  modelo<8,8,string> modelo1;
  modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
  modelo1.llenar();


  //parámetros
  modelo1.set_temp(0);
  modelo1.set_kb(0);
  modelo1.set_A_prob(1);
  modelo1.set_influencia_externa(10);
  modelo1.set_condicion_externa("omnivoro");
  modelo1.set_influencia_primeros_vecinos(influencia_vecinos);
  influencia_vecinos[0][0] = 10;
  influencia_vecinos[1][1] = 0;
  influencia_vecinos[2][2] = 20;
  const int longitud_paso=1;
  const int numero_pasos=100;

  ofstream energia_t,
           entropia_t,
	   proporciones_t;
  energia_t.open("energia_t.txt");
  entropia_t.open("entropia_t.txt");
  proporciones_t.open("proporciones_t.txt");

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
			proporciones_t<<num_estados[j]<<"\t";
		}
		proporciones_t<<"\n";
		energia_t<<i<<"\t"<<modelo1.energia()<<"\n";
                k=0;
	}  
  }

  energia_t.close();
  entropia_t.close();
  proporciones_t.close();

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
