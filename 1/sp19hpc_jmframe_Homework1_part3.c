#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

int main () {
      
   int l, m, n, i, j, k;
   for (l = 48; l <= 1536; l += l){ 
      n=l;
      // Number of total elements in the matrices, also the length of the 
      // one dimensional vectors to do the matix multiplications.
      m = n*n;
      printf("The matrix will be %d by %d, with a total of %d elements \n", n,n,m);
    
      //setting up an array with x rows and y columns
      double *A, *B, *C0, *C3;
      A = malloc(m * sizeof *A);
      B = malloc(m * sizeof *B);
      C0 = malloc(m * sizeof *C0); // Brute
      C3 = malloc(m * sizeof *C3); // dgemm3
      /////////////////////////////////////////////
      // SETTING UP THE MATRICES FOR MULTIPLICATION
      /////////////////////////////////////////////
      srand ( time ( NULL));
      for ( i = 0; i < n; i++ ) {
         for ( j = 0; j < n; j++ ){
            A[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
            B[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
            C0[i * n + j] = 0;
            C3[i * n + j] = 0;
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
      clock_t start_3 = clock();
  
      //////////////////////////////////////////////////////////////////
      // MAIN MATRIX MULTIPLICATION LOOP WITH OPTIMIZED REGISTER MEMORY
      /////////////////////////////////////////////////////////////////
      for ( i = 0; i < n; i += 3 ) 
         for ( j = 0; j < n; j += 3 ) {
            register double c00 = C3[i * n + j];
            register double c10 = C3[(i + 1) * n + j];
            register double c20 = C3[(i + 2) * n + j];
            register double c01 = C3[i * n + (j + 1)];
            register double c11 = C3[(i + 1) * n + (j + 1)];
            register double c21 = C3[(i + 2) * n + (j + 1)];
            register double c02 = C3[i * n + (j + 2)];
            register double c12 = C3[(i + 1) * n + (j + 2)];
            register double c22 = C3[(i + 2) * n + (j + 2)];
            for ( k = 0; k < n; k += 3 ) {
               register double a00 = A[i * n + k];
               register double a10 = A[(i + 1) * n + k];
               register double a20 = A[(i + 2) * n + k];
               register double b00 = B[k * n + j]; 
               register double b01 = B[k * n + (j + 1)]; 
               register double b02 = B[k * n + (j + 2)];
               c00 = c00 + a00 * b00;
               c10 = c10 + a10 * b00;
               c20 = c20 + a20 * b00;
               c01 = c01 + a00 * b01;
               c11 = c11 + a10 * b01;
               c21 = c21 + a20 * b01; 
               c02 = c02 + a00 * b02;
               c12 = c12 + a10 * b02;
               c22 = c22 + a20 * b02;
               // These variable keep the same name, but no longer represent these values.
               a00 = A[i * n + (k + 1)];
               a10 = A[(i + 1) * n + (k + 1)];
               a20 = A[(i + 2) * n + (k + 1)];
               b00 = B[(k + 1) * n + j];
               b01 = B[(k + 1) * n + (j + 1)];
               b02 = B[(k + 1) * n + (j + 2)];
               // Now add the values to the C matrix
               c00 = c00 + a00 * b00;
               c10 = c10 + a10 * b00;
               c20 = c20 + a20 * b00;
               c01 = c01 + a00 * b01;
               c11 = c11 + a10 * b01;
               c21 = c21 + a20 * b01; 
               c02 = c02 + a00 * b02;
               c12 = c12 + a10 * b02;
               c22 = c22 + a20 * b02;
               // Again, these values do not correspond with the variable names.
               a00 = A[i * n + (k + 2)];
               a10 = A[(i + 1) * n + (k + 2)];
               a20 = A[(i + 2) * n + (k + 2)];
               b00 = B[(k + 2) * n + j];
               b01 = B[(k + 2) * n + (j + 1)];
               b02 = B[(k + 2) * n + (j + 2)];
               // Make the final sum into the resulting C matrix
               c00 = c00 + a00 * b00;
               c10 = c10 + a10 * b00;
               c20 = c20 + a20 * b00;
               c01 = c01 + a00 * b01;
               c11 = c11 + a10 * b01;
               c21 = c21 + a20 * b01; 
               c02 = c02 + a00 * b02;
               c12 = c12 + a10 * b02;
               c22 = c22 + a20 * b02;
            }
            // Not move the data from the register to slow memory
            C3[i * n + j] = c00;
            C3[(i + 1) * n + j] = c10;
            C3[(i + 2) * n + j] = c20;
            C3[i * n + (j + 1)] = c01;
            C3[(i + 1) * n + (j + 1)] = c11;
            C3[(i + 2) * n + (j + 1)] = c21;
            C3[i * n + (j + 2)] = c02;
            C3[(i + 1) * n + (j + 2)] = c12;
            C3[(i + 2) * n + (j + 2)] = c22;
         }
      
      // recording the time of completion for the matrix multiplication
      clock_t end_3 = clock();
      double total_3 = ((double)(end_3 - start_3)) / CLOCKS_PER_SEC;
      printf("dgemm3\n");
      printf("Total matrix multiplication run time... \n");
      printf(" with MAXIMIZED register reuse: %f\n", total_3);   
      printf("Performance in GFLOPS: %f\n", (2*n^3)/total_3);   
      printf("\n");
      
      // Check that the matrix multiplactions are the same for each
      // NOTE: the matrix multiplications were tested in dgemm1 & dgemm0_1D
      double maxDiff3 = 0;
      double diff3;
      for ( i = 0; i < m; i++ ) {
         diff3 = C3[i] - C0[i];
         if ( fabs(diff3) > maxDiff3 ){
            maxDiff3 = diff3;
            printf("ELEMENTS IN MATRIX C0 & C3 ARE NOT THE SAME\t C0[%d] = %f\t C3[%d] = %f\n", i,C0[i], i, C3[i]);
         }     
      } 
      printf("The maximum difference between dgemm0 and dgemm3 is: %lf\n",maxDiff3);
      printf("-------------------------------------------------------------\n\n\n");
   }
   return 0;
}
//End of program 
