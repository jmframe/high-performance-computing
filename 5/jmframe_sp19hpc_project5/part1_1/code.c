#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#define N 2000

//Printing the Matrices to make sure the code worked.
void printMatrices(double **A, double **B, double **C, int i, int j, int nPrint){
   printf("Matrix A: \n");
   for (i = 0; i < nPrint; i++) {
     for (j = 0; j < nPrint; j++) {
       printf("%f\t", A[i][j]);
     }
     printf("\n");
   }
   printf("Matrix B: \n");
   for (i = 0; i < nPrint; i++) {
     for (j = 0; j < nPrint; j++) {
       printf("%f\t", B[i][j]);
     }
     printf("\n");
   }
   printf("Matrix C: \n");
   for (i = 0; i < nPrint; i++) {
     for (j = 0; j < nPrint; j++) {
       printf("%f\t", C[i][j]);
     }
     printf("\n");
   }
}


int main(int argc, char *argv){
   double sum, threadSum = 0;
   double start, end; //used for timing
   int threadz = 1;
   int threadCount = 0;
   int mCount = 0;
   int i, j, k;
   // Define the arrays as vectors with dynamic memory
   double **A, **B, **C;
   A = (double**)malloc(N * sizeof(double*));
   B = (double**)malloc(N * sizeof(double*));
   C = (double**)malloc(N * sizeof(double*));
   for(i = 0; i < N; i++){
      A[i] = (double*)malloc(N * sizeof(double));
      B[i] = (double*)malloc(N * sizeof(double));
      C[i] = (double*)malloc(N * sizeof(double));
   }
   
   //Assigning values to the matrices.
   //This is an easy way to assess the correctness of the calculation
   //Since the A and B matrices are set to 1.000, all results in the 
   //C matrix should be equal to the value of N
   for (i = 0; i < N; i++){
      for (j = 0; j < N; j++){
         A[i][j] = 1;
         B[i][j] = 1;
      }
   }
   while(threadz <= 32){

      omp_set_num_threads(threadz); //set number of threads here

      //Re-setting values to the C matrix.
      for (i = 0; i < N; i++){
         for (j = 0; j < N; j++){
            C[i][j] = 0;
         }
      }
   
      //Main multiplication loop
      start = omp_get_wtime(); //start time measurement
   
      ///////////////////////// SET OPEN MP FOR LOOP //////////////////
      #pragma omp parallel for shared(A,B,C) private(i, j, k, sum)
      for (i = 0; i < N; i++){
         for (j = 0; j < N; j++){
   
            sum=0;
            for (k = 0; k < N; k++){
               sum += A[i][k] * B[k][j];
            }
   
            // Assignming sum to C matrix
            C[i][j] = sum;

         }
      }
      
      end = omp_get_wtime(); //end time measurement
      threadSum = threadSum + (end-start);
      printf("For the number of threads %d, calculation %d of 5\n", threadz, threadCount+1);
      printf("Time of computation: %f seconds\n", end - start);
      
      // Print to debug
      int nPrint = 2; //print the firsr # rows/columns of the matrices.
      printMatrices(A, B, C, i, j, nPrint);

      // testing the matrix multiplication for accuracy
      // The A & B matrices are set to 1
      // so the multiplication matrix (C) should have values of N everywhere.
      for(i = 0; i < N; i++){
         for(j = 0; j < N; j++){
            if(C[i][j] != N){
               printf("There is an error in the matrix-matrix multiplication!"); 
               printf("The element C[%d][%d] = %f\n", i, j, C[i][j]);
               exit(1);
            }
         }
      }

      // Now to get the average results, move the thread count forward after 5 times
      threadCount = threadCount + 1;  
      if (threadCount > 4){
         printf("The average time of computation for %d threads is: %f seconds\n\n\n", threadz, threadSum/threadCount);
         printf("-------------------------------------------------------------\n\n");
         threadz = threadz * 2;
         threadCount = 0;
         threadSum = 0;
      }

   }

   printf("Program finished successfully");
   return 0;

}
