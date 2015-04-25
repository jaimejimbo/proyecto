/*
 *
 *
 *
 */

#include <iostream>
using namespace std;
#include "metropolis.hpp"


int main(int argc, char** argv){

  modelo<10,10> modelo1;
  const int N_posibles_estados = 2;  
  double posibles_estados[N_posibles_estados];
  posibles_estados[0] = -1/2;
  posibles_estados[1] = 1/2;
  modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
  modelo1.llenar();
  cout << modelo1(1,2);
  return 0;
}



