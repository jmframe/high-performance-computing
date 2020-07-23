#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define N 1040
#define T 100000

//Printing the Matrices to make sure the code worked.
void printMatrix(double ***M, int i, int j, int t){
   int nprint = N/8;
   printf("The result af %d time step: \n",t);
   for (i = 0; i < N; i += nprint) {
      for (j = 0; j < N; j += nprint) {
         printf("%f\t", M[i][j][0]);
      }
      printf("%f\t", M[i][N-1][0]);
      printf("\n");
   }
   for (j = 0; j < N; j += nprint) {
      printf("%f\t", M[N-1][j][0]);
   } 
   printf("%f\t", M[N-1][N-1][0]);
   printf("\n");
}

//compare single thread to make sure multi-thread code worked.
int checkCalc(double ***M, int i, int j){
   float single, multi;
   int diff = 0, maxDiff = 0;
   for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
         single = round(M[i][j][2]);
         multi = round(M[i][j][0]);
         diff = abs(single - multi);
         if(abs(single-multi) > diff){
            maxDiff = abs(single - multi);
         }
      }
   }
   return maxDiff;
}

int main(int argc, char *argv){
   double start, end; //used for timing
   double runTime, seqRunTime, speedup; //performance analysis values
   int threadz = 1; //number of threads for parallel loop
   int d = 2;       //the third demension of the array
   int i, j, k, t;  //just counting variables
   int a, b, c; //variables to define demention location
   int numPrints = 5; //how many times to print the matrix
   int divisor = T/numPrints;// the divisor to print the appropriate matrices
   int error = 0; // a value to assess the error in the multi-thread calc.

   // Define the arrays as vectors with dynamic memory
   double *** h = (double ***)malloc(N * sizeof(double**));
   for (i = 0; i< N; i++) {
      h[i] = (double **) malloc(N * sizeof(double *));
      for (j = 0; j < N; j++) {
         h[i][j] = (double *)malloc(d * sizeof(double));
      }
   }

   // Loop through the calculation for different thread values
   for(threadz = 1; threadz < 33; threadz = threadz * 2){ 

      //Assigning values to the matrices.
      for (i = 0; i < N; i++){
         for (j = 0; j < N; j++){
            for (k = 0; k < d; k++){
               h[i][j][k] = 20;
               if(i == 0 && N * 0.3 <= j && j <= N * 0.7){
                  h[i][j][k] = 100;
               }
            }
         }
      }

      //setting the matrix location variables
      a=0;
      b=1;

      start = omp_get_wtime(); //start time measurement
   
      // MAIN CALCULATION LOOP ////////////////////////////
      for (t = 0; t < T; t++) { 
         //   Set the Open MP stuff
         omp_set_num_threads(threadz); //set number of threads here
         ///////////////////////// SET OPEN MP FOR LOOP //////////////////
         #pragma omp parallel for shared(h) private(i, j)
         for (i = 1; i < (N-1); i++){
            for (j = 1; j < (N-1); j++){
               h[i][j][a] = (h[i-1][j][b] + h[i+1][j][b] + h[i][j-1][b] + h[i][j+1][b]) / 4; 
            }
         }
         c = a;
         a = b;
         b = c;
         
         // Print every so often to see the progress of the heat distribution 
         if(t % divisor == 0){
            printf("h matrix at time: %d\n",t);
            printMatrix(h, i, j, t);
         }

      } //END THE MAIN CALCULATION LOOP
   
      end = omp_get_wtime(); //end time measurement
      runTime = end - start;

      //Check that the multi-thread calc is same as single.
      if(threadz == 1){
         seqRunTime = runTime;
         for (i = 0; i < N; i++){
            for (j = 0; j < N; j++){
               h[i][j][2] = h[i][j][0];
            }
         }
         printf("Single thread result\n");
      }
      if(threadz > 1){
         error = checkCalc(h, i, j);
         if(error >= 1){
            printf("Error for the number of threads: %d\n", threadz);
            printf("The max difference between this calc and single thread is: %d\n", error);
         }
         //calculate speedup while we are here.
         speedup = seqRunTime/(runTime/threadz);
      }
 
      printf("The final result after %d time steps: \n",T);
      printMatrix(h,i,j,T);  
      printf("For the number of threads %d\n", threadz);
      printf("The total run time is: %f\n", runTime );
      if(threadz > 1){
        printf("The speedup of this calculation is: %f\n", speedup);
      }
      printf("------------------------------------------------------------------\n\n\n");
      error = 0;

   }

   printf("Program finished successfully");
   return 0;

}
