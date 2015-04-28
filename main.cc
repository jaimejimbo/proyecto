/*
 *
 *
 *
 */

#include <iostream>
using namespace std;
#include "metropolis.hpp"



int main(int argc, char** argv){
  const int N_posibles_estados = 3;
  string posibles_estados[N_posibles_estados];
  posibles_estados[0] = "omnivoro";
  posibles_estados[1] = "carnivoro";
  posibles_estados[2] = "vegetariano";
  double *influencia_vecinos = new double[3];
  influencia_vecinos[0] = 10000;
  influencia_vecinos[1] = 0;
  influencia_vecinos[2] = 0;
  modelo<8,8,string> modelo1;
  modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
  modelo1.llenar();
  modelo1.set_temp(0);
  modelo1.set_kb(0);
  modelo1.set_A_prob(1);
  modelo1.set_influencia_externa(10000);
  modelo1.set_condicion_externa("omnivoro");
  modelo1.set_influencia_primeros_vecinos(influencia_vecinos);
  cout<<modelo1.energia()<<"\n";
  for (int i=0; i<100; i++) modelo1.cambiar_estado();
  cout<<modelo1.energia()<<"\n";
  cout<<modelo1.to_string();
}



