//////////////////////     INCLUDE LIBRARIES    ////////////////////
#include <stdio.h>

//Only use one of these
#include <stdlib.h>
//#include <string.h>

#include <stdbool.h>
#include <time.h>
////////////////////

////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void dgemm0(double *A, double *B, double *C0, int n){
   int i, j, k;
   for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < n; j++ ) {
         for ( k = 0; k < n; k++ ) {
           C0[i * n + j] += A[i * n + k] * B[k * n + j];
         }
      }
   }
}

////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void ijk_blocked(double *A, double *B, double *C1, int n, int Blc){
   int i, j, k, i1, j1, k1;
   for ( i = 0; i < n; i += Blc ) {
      for ( j = 0; j < n; j += Blc ) {
         for ( k = 0; k < n; k += Blc ) { 
            // B x B mini matrix multiplications
            for ( i1 = i; i1 < (i + Blc); i1 += 2 ) { 
               for ( j1 = j; j1 < (j + Blc); j1 += 2 ) {
                  register double C00 = C1[(i1 + 0) * n + (j1 + 0)];
                  register double C01 = C1[(i1 + 0) * n + (j1 + 1)];
                  register double C10 = C1[(i1 + 1) * n + (j1 + 0)];
                  register double C11 = C1[(i1 + 1) * n + (j1 + 1)];
                  for ( k1 = k; k1 < (k + Blc); k1 += 2 ) {
                     register double A00 = A[(i1 + 0) * n + (k1 + 0)];
                     register double A01 = A[(i1 + 0) * n + (k1 + 1)];
                     register double A10 = A[(i1 + 1) * n + (k1 + 0)];
                     register double A11 = A[(i1 + 1) * n + (k1 + 1)];
                     register double B00 = B[(k1 + 0) * n + (j1 + 0)];
                     register double B01 = B[(k1 + 0) * n + (j1 + 1)];
                     register double B10 = B[(k1 + 1) * n + (j1 + 0)];
                     register double B11 = B[(k1 + 1) * n + (j1 + 1)];
                     C00 += A00 * B00 + A01 * B10;
                     C01 += A00 * B01 + A01 * B11;
                     C10 += A10 * B00 + A11 * B10;
                     C11 += A10 * B01 + A11 * B11;
                  }
                  C1[(i1 + 0) * n + (j1 + 0)] = C00;
                  C1[(i1 + 0) * n + (j1 + 1)] = C01;
                  C1[(i1 + 1) * n + (j1 + 0)] = C10;
                  C1[(i1 + 1) * n + (j1 + 1)] = C11;
               }
            }
         }
      }
   }
}
////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void jik_blocked(double *A, double *B, double *C1, int n, int Blc){
   int i, j, k, i1, j1, k1;
   for ( j = 0; j < n; j += Blc ) {
      for ( i = 0; i < n; i += Blc ) {
         for ( k = 0; k < n; k += Blc ) { 
            // B x B mini matrix multiplications
            for ( j1 = j; j1 < (j + Blc); j1 += 2 ) { 
               for ( i1 = i; i1 < (i + Blc); i1 += 2 ) {
                  register double C00 = C1[(i1 + 0) * n + (j1 + 0)];
                  register double C01 = C1[(i1 + 0) * n + (j1 + 1)];
                  register double C10 = C1[(i1 + 1) * n + (j1 + 0)];
                  register double C11 = C1[(i1 + 1) * n + (j1 + 1)];
                  for ( k1 = k; k1 < (k + Blc); k1 += 2 ) {
                     register double A00 = A[(i1 + 0) * n + (k1 + 0)];
                     register double A01 = A[(i1 + 0) * n + (k1 + 1)];
                     register double A10 = A[(i1 + 1) * n + (k1 + 0)];
                     register double A11 = A[(i1 + 1) * n + (k1 + 1)];
                     register double B00 = B[(k1 + 0) * n + (j1 + 0)];
                     register double B01 = B[(k1 + 0) * n + (j1 + 1)];
                     register double B10 = B[(k1 + 1) * n + (j1 + 0)];
                     register double B11 = B[(k1 + 1) * n + (j1 + 1)];
                     C00 += A00 * B00 + A01 * B10;
                     C01 += A00 * B01 + A01 * B11;
                     C10 += A10 * B00 + A11 * B10;
                     C11 += A10 * B01 + A11 * B11;
                  }
                  C1[(i1 + 0) * n + (j1 + 0)] = C00;
                  C1[(i1 + 0) * n + (j1 + 1)] = C01;
                  C1[(i1 + 1) * n + (j1 + 0)] = C10;
                  C1[(i1 + 1) * n + (j1 + 1)] = C11;
               }
            }
         }
      }
   }
}
////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void kij_blocked(double *A, double *B, double *C1, int n, int Blc){
   int i, j, k, i1, j1, k1;
   for ( k = 0; k < n; k += Blc ) {
      for ( i = 0; i < n; i += Blc ) {
         for ( j = 0; j < n; j += Blc ) { 
            // B x B mini matrix multiplications
            for ( k1 = k; k1 < (k + Blc); k1 += 2 ) { 
               for ( i1 = i; i1 < (i + Blc); i1 += 2 ) {
                  register double A00 = A[(i1 + 0) * n + (k1 + 0)];
                  register double A01 = A[(i1 + 0) * n + (k1 + 1)];
                  register double A10 = A[(i1 + 1) * n + (k1 + 0)];
                  register double A11 = A[(i1 + 1) * n + (k1 + 1)];
                  for ( j1 = j; j1 < (j + Blc); j1 += 2 ) {
                     register double B00 = B[(k1 + 0) * n + (j1 + 0)];
                     register double B01 = B[(k1 + 0) * n + (j1 + 1)];
                     register double B10 = B[(k1 + 1) * n + (j1 + 0)];
                     register double B11 = B[(k1 + 1) * n + (j1 + 1)];
                     C1[(i1 + 0) * n + (j1 + 0)] += A00 * B00 + A01 * B10;
                     C1[(i1 + 0) * n + (j1 + 1)] += A00 * B01 + A01 * B11;
                     C1[(i1 + 1) * n + (j1 + 0)] += A10 * B00 + A11 * B10;
                     C1[(i1 + 1) * n + (j1 + 1)] += A10 * B01 + A11 * B11;
                  }
               }
            }
         }
      }
   }
}
////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void ikj_blocked(double *A, double *B, double *C1, int n, int Blc){
   int i, j, k, i1, j1, k1;
   for ( i = 0; i < n; i += Blc ) {
      for ( k = 0; k < n; k += Blc ) {
         for ( j = 0; j < n; j += Blc ) { 
            // B x B mini matrix multiplications
            for ( i1 = i; i1 < (i + Blc); i1 +=2 ) { 
               for ( k1 = k; k1 < (k + Blc); k1 += 2 ) {
                  register double A00 = A[(i1 + 0) * n + (k1 + 0)];
                  register double A01 = A[(i1 + 0) * n + (k1 + 1)];
                  register double A10 = A[(i1 + 1) * n + (k1 + 0)];
                  register double A11 = A[(i1 + 1) * n + (k1 + 1)];
                  for ( j1 = j; j1 < (j + Blc); j1 +=2 ) {
                     register double B00 = B[(k1 + 0) * n + (j1 + 0)];
                     register double B01 = B[(k1 + 0) * n + (j1 + 1)];
                     register double B10 = B[(k1 + 1) * n + (j1 + 0)];
                     register double B11 = B[(k1 + 1) * n + (j1 + 1)];
                     C1[(i1 + 0) * n + (j1 + 0)] += A00 * B00 + A01 * B10;
                     C1[(i1 + 0) * n + (j1 + 1)] += A00 * B01 + A01 * B11;
                     C1[(i1 + 1) * n + (j1 + 0)] += A10 * B00 + A11 * B10;
                     C1[(i1 + 1) * n + (j1 + 1)] += A10 * B01 + A11 * B11;
                  }
               }
            }
         }
      }
   }
}
////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void jki_blocked(double *A, double *B, double *C1, int n, int Blc){
   int i, j, k, i1, j1, k1;
   for ( j = 0; j < n; j += Blc ) {
      for ( k = 0; k < n; k += Blc ) {
         for ( i = 0; i < n; i += Blc ) { 
            // B x B mini matrix multiplications
            for ( j1 = j; j1 < (j + Blc); j1 += 2 ) { 
               for ( k1 = k; k1 < (k + Blc); k1 += 2 ) {
                  register double B00 = B[(k1 + 0) * n + (j1 + 0)];
                  register double B01 = B[(k1 + 0) * n + (j1 + 1)];
                  register double B10 = B[(k1 + 1) * n + (j1 + 0)];
                  register double B11 = B[(k1 + 1) * n + (j1 + 1)];
                  for ( i1 = i; i1 < (i + Blc); i1 += 2 ) {
                     register double A00 = A[(i1 + 0) * n + (k1 + 0)];
                     register double A01 = A[(i1 + 0) * n + (k1 + 1)];
                     register double A10 = A[(i1 + 1) * n + (k1 + 0)];
                     register double A11 = A[(i1 + 1) * n + (k1 + 1)];
                     C1[(i1 + 0) * n + (j1 + 0)] += A00 * B00 + A01 * B10;
                     C1[(i1 + 0) * n + (j1 + 1)] += A00 * B01 + A01 * B11;
                     C1[(i1 + 1) * n + (j1 + 0)] += A10 * B00 + A11 * B10;
                     C1[(i1 + 1) * n + (j1 + 1)] += A10 * B01 + A11 * B11;
                  }
               }
            }
         }
      }
   }
}
////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void kji_blocked(double *A, double *B, double *C1, int n, int Blc){
   int i, j, k, i1, j1, k1;
   for ( k = 0; k < n; k += Blc ) {
      for ( j = 0; j < n; j += Blc ) {
         for ( i = 0; i < n; i += Blc ) { 
            // B x B mini matrix multiplications
            for ( k1 = k; k1 < (k + Blc); k1 += 2 ) { 
               for ( j1 = j; j1 < (j + Blc); j1 += 2 ) {
                  register double B00 = B[(k1 + 0) * n + (j1 + 0)];
                  register double B01 = B[(k1 + 0) * n + (j1 + 1)];
                  register double B10 = B[(k1 + 1) * n + (j1 + 0)];
                  register double B11 = B[(k1 + 1) * n + (j1 + 1)];
                  for ( i1 = i; i1 < (i + Blc); i1 += 2 ) {
                     register double A00 = A[(i1 + 0) * n + (k1 + 0)];
                     register double A01 = A[(i1 + 0) * n + (k1 + 1)];
                     register double A10 = A[(i1 + 1) * n + (k1 + 0)];
                     register double A11 = A[(i1 + 1) * n + (k1 + 1)];
                     C1[(i1 + 0) * n + (j1 + 0)] += A00 * B00 + A01 * B10;
                     C1[(i1 + 0) * n + (j1 + 1)] += A00 * B01 + A01 * B11;
                     C1[(i1 + 1) * n + (j1 + 0)] += A10 * B00 + A11 * B10;
                     C1[(i1 + 1) * n + (j1 + 1)] += A10 * B01 + A11 * B11;
                  }
               }
            }
         }
      }
   }
}

////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void checkresults(double *C0, double *C1, int m, clock_t start_t, clock_t end_t){
   int i;
   // Check that the matrix multiplactions are the same for each
   // NOTE: the matrix multiplications were tested in dgemm0 & Blocked
   double maxDiff = 0;
   double diff;
   for ( i = 0; i < m; i++ ) {
      diff = C1[i] - C0[i];
      if ( abs(diff) > maxDiff ){
         maxDiff = diff;
         printf("ELEMENTS IN MATRIX C0 & C1 ARE NOT THE SAME\t C0[%d] = %f\t C1[%d] = %f\n", i,C0[i], i, C1[i]);
      }
   }
   double total_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
   printf("Total matrix multiplication run time: %f seconds\n", total_t); 
   printf("The maximum difference between blocked and unblocked is: %lf\n",maxDiff);
   printf("-------------------------------------------------------------\n");
}
////////////////////////////////////   VOID FUNCTION   ////   VOID FUNCTION   ////////////////////////////////
void clearmatrixC(double *C1, int n){ 
   int i, j;
   // Clear the C matrices, to make sure we dont cross over results.
   for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < n; j++ ){
         C1[i * n + j] = 0;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN Program, which sets up and executes ////////////////////////////////////////////////////////////////
// matrix multibplications for blocked and /////////////////////////////////////////////////////////////////
// unblocked algorithms.////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main () {
   int m, n, i, j, k, Blc;
   clock_t start_t, end_t;

   // MATRIX SIZE
   n=512;
   // Number of total elements in the matrices, also the length of the 
   // one dimensional vectors to do the matix multiplications.
   m = n*n;
   printf("----------------------------------------\n");
   printf("The matrix size is %d by %d \n", n,n);
   
   //setting up an array with x rows and y columns
   double *A, *B, *C0, *C1;
   A = malloc(m * sizeof *A);
   B = malloc(m * sizeof *B);
   C0 = malloc(m * sizeof *C0);
   C1 = malloc(m * sizeof *C1); 

   /////////////////////////////////////////////
   // SETTING UP THE MATRICES FOR MULTIPLICATION
   for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < n; j++ ){
         A[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
         B[i * n + j] = (double)rand()/RAND_MAX*2.0-1.0;
         C0[i * n + j] = 0;
         C1[i * n + j] = 0;
      }
   }
 
   /////////////////////////////////////////////
   // RUN THE BASIC MULTIPLICATION SCRIPT
   // TO ENSURE RESULTS OF BLOCK MULTIPLICATION
   // Time this program, to use for efficiency tests
   dgemm0(A, B, C0, n);

   printf("--------------------------------------------------\n");
   printf("Looping through block sizes.\n");
   printf("Running each algorithm for each block size\n");
   printf("--------------------------------------------------\n\n\n\n\n");

   ///////////////////////////////////////////////////////////////////////
   // LOOP THROUGH THE BLOCK SIZE TO FIND OPTIMIZE
   //////////////////////////////////////////////////////////////////////
   for ( Blc = 8; Blc <=128; Blc += Blc ) {
      ///////////////////////////////////////////////////////////////////
      // One of these for each algorithm
      printf("BLOCK SIZE: %d, ALGORITHM: ijk\n", Blc);
      start_t = clock();
      ijk_blocked(A, B, C1, n, Blc);
      end_t = clock();
      checkresults(C0, C1, m, start_t, end_t);
      clearmatrixC(C1, n);
      printf("BLOCK SIZE: %d, ALGORITHM: ikj\n", Blc);
      start_t = clock();
      ikj_blocked(A, B, C1, n, Blc);
      end_t = clock();
      checkresults(C0, C1, m, start_t, end_t);
      clearmatrixC(C1, n);
      printf("BLOCK SIZE: %d, ALGORITHM: jik\n", Blc);
      start_t = clock();
      jik_blocked(A, B, C1, n, Blc);
      end_t = clock();
      checkresults(C0, C1, m, start_t, end_t);
      clearmatrixC(C1, n);
      printf("BLOCK SIZE: %d, ALGORITHM: jki\n", Blc);
      start_t = clock();
      jki_blocked(A, B, C1, n, Blc);
      end_t = clock();
      checkresults(C0, C1, m, start_t, end_t);
      clearmatrixC(C1, n);
      printf("BLOCK SIZE: %d, ALGORITHM: kij\n", Blc);
      start_t = clock();
      kij_blocked(A, B, C1, n, Blc);
      end_t = clock();
      checkresults(C0, C1, m, start_t, end_t);
      clearmatrixC(C1, n);
      printf("BLOCK SIZE: %d, ALGORITHM: kji\n", Blc);
      start_t = clock();
      kji_blocked(A, B, C1, n, Blc);
      end_t = clock();
      checkresults(C0, C1, m, start_t, end_t);
      clearmatrixC(C1, n);
      printf("---------------------------------------------------------------------\n\n");
   }
   printf("--- END PROGRAM -------------------------------------------------\n");
   printf("------------------END PROGRAM -----------------------------------\n");
   printf("-------------------------------END PROGRAM-----------------------");
   return 0;
}
//End of program 
