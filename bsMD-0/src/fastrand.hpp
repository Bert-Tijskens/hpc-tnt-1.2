static int g_seed = 0;

inline int fastrand() { 
  g_seed = (214013*g_seed+2531011); 
  return (g_seed>>16)&0x7FFF; 
} 

#include <limits>

template <class T>
T fast_uniform_real() {
    int i = fastrand();
    T tmax = T( std::numeric_limits<int>::max() );
    T ti   = T(i);
    T t    = ti/tmax;
    return t;    
}