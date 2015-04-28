/*
 *
 *
 *
 */

#include <iostream>
using namespace std;
#include "metropolis.hpp"



int main(int argc, char** argv){
  string posibles_estados[2];
  int N_posibles_estados = 2;
  posibles_estados[0] = "+";
  posibles_estados[1] = "-";
  modelo<10,10,string> modelo1;
  /*modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
  modelo1.llenar();
  modelo1.set_temp(0);
  modelo1.set_kb(0);
  modelo1.set_A_prob(1);
  modelo1.set_influencia_externa(0);
  modelo1.set_condicion_externa(0);
  modelo1.set_influencia_primeros_vecinos(1);
  cout<<modelo1.energia()<<"\n";
  modelo1.cambiar_estado();
  cout<<modelo1.energia()<<"\n";*/

      
}



