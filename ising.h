
class Matriz //de orden
{
  public:
    Matriz(unsigned int N_filas, unsigned int N_columnas);
  private:
    double** data;
    unsigned int N_filas,
                 N_columnas;
};

Matriz::Matriz(unsigned int f, unsigned int c)
{
  data = new double[f][c]
  N_filas = f;
  N_columnas = c;
}