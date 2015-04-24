
#include "ising.h"

Matriz::Matriz(const unsigned int f, const unsigned int c)
{
  data = new double[f][c]
  N_filas = f;
  N_columnas = c;
}
