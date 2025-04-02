#include <stdio.h> 
#include <omp.h>

void main ()  { 
    double start;
    double end;
     start = omp_get_wtime();
	     for (int i=0;i<100000;i++) {} 
     end = omp_get_wtime();
     printf("Praca zajela %f sekund\n", end - start);

     // omp_get_wtick()  - precyzja pomiaru czasu
     printf("Precyzja pomiaru: %.9f\n",omp_get_wtick());
     // Precyzja pomiaru: 0.000000446[s]
}