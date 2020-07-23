#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

int main () {
      
   int l, m, n, i, j, k;

   // LOOP THROUGH TEST VALUES
   for ( l = 8; l <= 2048; l += l ){
      n=l;
      // Number of total elements in the matrices, also the length of the 
      // one dimensional vectors to do the matix multiplications.
      m = n*n;
      printf("The matrix will be %d by %d, with a total of %d elements \n", n,n,m);
    
      //setting up an array with x rows and y columns
      double *A, *B, *C0, *C1, *C2;
      A = malloc(m * sizeof *A);
      B = malloc(m * sizeof *B);
      C0 = malloc(m * sizeof *C0); // Brute
      C1 = malloc(m * sizeof *C1); // dgemm1
      C2 = malloc(m * sizeof *C2); // dgemm2
      /////////////////////////////////////////////
      // SETTING UP THE MATRICES FOR MULTIPLICATION
      /////////////////////////////////////////////
      for ( i = 0; i < n; i++ ) {
         for ( j = 0; j < n; j++ ){
            A[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
            B[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
            C0[i * n + j] = 0;
            C1[i * n + j] = 0;
            C2[i * n + j] = 0;
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
      printf("Performance in GFLOPS: %f\n", (2*n^3)/total_0);   
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
      printf("Performance in GFLOPS: %f\n", (2*n^3)/total_1);   
      printf("\n");
    
      // Time this program, to use for efficiency tests
      clock_t start_2 = clock();
  
      //////////////////////////////////////////////////////////////////
      // MAIN MATRIX MULTIPLICATION LOOP WITH AGGRESIVE REGISTER MEMORY
      /////////////////////////////////////////////////////////////////
      for ( i = 0; i < n; i += 2 ) {
         for ( j = 0; j < n; j += 2 ) {
           register double c0 = C2[i * n + j];
           register double c1 = C2[(i + 1) * n + j];
           register double c2 = C2[i * n + (j + 1)];
           register double c3 = C2[(i + 1) * n + (j + 1)];
           for ( k = 0; k < n; k += 2 ) {
              register double a0 = A[i * n + k];
              register double a1 = A[i * n + (k + 1)];
              register double a2 = A[(i + 1) * n + k];
              register double a3 = A[(i + 1) * n + (k + 1)];
              register double b0 = B[k * n + j]; 
              register double b1 = B[(k + 1) * n + j];
              register double b2 = B[k * n + (j + 1)]; 
              register double b3 = B[(k + 1) * n + (j + 1)];
              c0 = a0 * b0 + a1 * b1 + c0;
              c1 = a2 * b0 + a3 * b1 + c1;
              c2 = a0 * b2 + a1 * b3 + c2;
              c3 = a2 * b2 + a3 * b3 + c3;
           }
           C2[i * n + j] = c0;
           C2[(i + 1) * n + j] = c1;
           C2[i * n + (j + 1)] = c2;
           C2[(i + 1) * n + (j + 1)] = c3;
         }
      }
 
      // recording the time of completion for the matrix multiplication
      clock_t end_2 = clock();
      double total_2 = ((double)(end_2 - start_2)) / CLOCKS_PER_SEC;
      printf("dgemm2\n");
      printf("Total matrix multiplication run time... \n");
      printf(" with AGGRESSIVE register reuse: %f\n", total_2);   
      printf("Performance in GFLOPS: %f\n", (2*n^3)/total_2);   
      printf("\n");

      // Check that the matrix multiplactions are the same for each
      // NOTE: the matrix multiplications were tested in dgemm1 & dgemm0_1D
      double maxDiff1 = 0;
      double maxDiff2 = 0;
      double diff1;
      double diff2;
      for ( i = 0; i < m; i++ ) {
         diff1 = C1[i] - C0[i];
         diff2 = C2[i] - C0[i];
         if ( abs(diff1) > maxDiff1 ){
            maxDiff1 = diff1;
         }
         if ( abs(diff2) > maxDiff2 ){
            maxDiff2 = diff2;
         }
      } 
      printf("The maximum difference between dgemm0 and dgemm1 is: %lf\n",maxDiff1);
      printf("The maximum difference between dgemm0 and dgemm2 is: %lf\n",maxDiff2);
      printf("-------------------------------------------------------------\n\n\n");
   }
   return 0;
}
//End of program 
