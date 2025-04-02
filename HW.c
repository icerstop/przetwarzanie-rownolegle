#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
int main(){
 #pragma omp parallel
  { 
   #ifdef _OPENMP
   printf("Hello, world. (tid=%d)\n",omp_get_thread_num()); 
   #else
   printf("No OpenMP\n");
   #endif 
  }
 return EXIT_SUCCESS;
}