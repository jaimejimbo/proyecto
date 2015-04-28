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
  posibles_estados[0] = "+";
  posibles_estados[1] = "0";
  posibles_estados[2] = "-";
  modelo<10,10,string> modelo1;
  modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
  modelo1.llenar();
  modelo1.set_temp(100);
  modelo1.set_kb(1);
  modelo1.set_A_prob(1);
  modelo1.set_influencia_externa(100);
  modelo1.set_condicion_externa("+");
  modelo1.set_influencia_primeros_vecinos(0);
  cout<<modelo1.energia()<<"\n";
  modelo1.cambiar_estado();
  cout<<modelo1.energia()<<"\n";
}



