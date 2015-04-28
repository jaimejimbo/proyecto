/*
 *
 *
 *
 */

#include <iostream>
using namespace std;
#include "metropolis.hpp"

typedef enum {estado1 = 1, estado2 = 2, estado3 = 3} posibles_estados;

int main(int argc, char** argv){
	modelo<10,10,posibles_estados> modelo1;
        int N_posibles_estados = 3;
        posibles_estados pe;
	modelo1.definir_posibles_estados(N_posibles_estados, pe);
	modelo1.llenar();
	modelo1.set_temp(0);
	modelo1.set_kb(0);
	modelo1.set_A_prob(1);
	modelo1.set_influencia_externa(0);
	modelo1.set_condicion_externa(0);
	modelo1.set_influencia_primeros_vecinos(1);
  cout<<modelo1.energia()<<"\n";
	modelo1.cambiar_estado();
  cout<<modelo1.energia()<<"\n";
        
}



