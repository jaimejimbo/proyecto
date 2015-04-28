/*
 *
 *
 *
 */



#define TEMPLATE                  template< size_t M, size_t N, typename T>



template< size_t N, typename T = float>
class estados{

public:
  estados();
  T get_estado(int indice);
private:
  T *posibles_estados = new T[N];

}


TEMPLATE
estados::estados()
{
  
}
