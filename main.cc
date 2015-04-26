/*
 *
 *
 *
 */

#include <iostream>
using namespace std;
#include "metropolis.hpp"


int main(int argc, char** argv){
	modelo<3,3> modelo1;
	int N_posibles_estados = 2;  
	double posibles_estados[N_posibles_estados];
	posibles_estados[0] = -1.0/2;
	posibles_estados[1] = 1.0/2;
	modelo1.definir_posibles_estados(N_posibles_estados, posibles_estados);
	modelo1.llenar();
	cout<<"Estado inicial:\n";
	modelo1.print();
	modelo1.cambiar_estado();
	cout<<"Estado final:\n";
  modelo1.print();
}



