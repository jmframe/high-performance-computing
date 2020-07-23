#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

int main () {
   
   //setting up an array with x rows and y columns
   int l, m, n, i ,j, k;
   for ( l = 64; l <= 2048; l += l ){
      //printf("Enter size of matrix:");
      //scanf("%d", &n);
      n=l;
      m = n * n;
      printf("The matrix will be %d by %d, with a total of %d elements \n", n,n,m);
     
      double random_value;
      srand ( time ( NULL));
    
      // Set up the vectors with random values.
      double A[m], B[m], C0[m], C1[m];
      for ( i = 0; i < n; i++ ) {
         for ( j = 0; j < n; j++ ){
            A[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
            B[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
            C0[i * n + j] = 0;
         }
      }
   
      // Time this program, to use for efficiency tests
      clock_t start_0 = clock(); 
   
      /////////////////////////////////////////////////////////////////
      // MAIN MATRIX MULTIPLICATION LOOP WITHOUT USING REGISTER MEMORY
      /////////////////////////////////////////////////////////////////
      for ( i = 0; i < n; i++ ) {
         for ( j = 0; j < n; j++ ) {
            for ( k = 0; k < n; k++ ) {
               C0[i * n + j] += A[i * n + k] * B[k * n+ j];
           }
         }
      }
   
      // recording the time of completion for the matrix multiplication
      clock_t end_0 = clock();
      double total_0 = ((double)(end_0 - start_0)) / CLOCKS_PER_SEC;
      printf("dgemm0\n");
      printf("Total matrix multiplication run time WITHOUT Register memory: %f\n", total_0);   
      printf("\n");
   
   
      // Time this program, to use for efficiency tests
      clock_t start_1 = clock();
   
      //////////////////////////////////////////////////////////////////
      // MAIN MATRIX MULTIPLICATION LOOP WITH REGISTER MEMORY
      /////////////////////////////////////////////////////////////////
      for ( i = 0; i < n; i++ ) {
         for ( j = 0; j < n; j++ ) {
           register double r = C1[i * n + j];
           for ( k = 0; k < n; k++ ) {
               r += A[i * n + k] * B[k * n + j];
           }
           C1[ i * n + j] = r;
         }
      }
   
      // recording the time of completion for the matrix multiplication
      clock_t end_1 = clock();
      double total_1 = ((double)(end_1 - start_1)) / CLOCKS_PER_SEC;
      printf("dgemm1\n");
      printf("Total matrix multiplication run time WITH Register memory: %f\n", total_1);   
      printf("\n");
    
      // Check that the matrix multiplactions are the same for each
      // NOTE: the matrix multiplications were tested in dgemm1 & dgemm0_1D
      double maxDiff = 0;
      double diff;
      for ( i = 0; i < m; i++ ) {
         diff = C1[m] - C0[m];
         if ( abs(diff) > maxDiff ){
            maxDiff = diff;
         }
      } 
      printf("The maximum difference between the two resulting matrices is: %lf\n",maxDiff);
      printf("-------------------------------------------------------------\n\n\n");
   }   
   return 0;
}
//End of program 
