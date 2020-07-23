// CS 491/591: High Performance Computing Project 3
// Parallel Sieve of Eratosthenes for Finding All Prime Numbers within 10^10

// Part 1 Simply modify the parallel Sieve of Eratosthenes program in class
// so that the program can find all prime members within 10^10 
// (hint: 10^10 is larger than value range of int type).

// Sieve of Eratosthenes 
// Programmed by Michael J. Quinn: 7 September 2001
// Code adapted by sp19hpc_jmframe (Jonathan Frame)
// March 2019

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   long long int  count;        // Local prime count
   double         elapsed_time; // Parallel execution time
   long long int  first;        // Index of first multiple
   long long int  global_count; // Global prime count
   long long int  high_value;   // Highest value on this proc
   long long int  i;
   int            id;           // Process ID number
   long long int  index;        // Index of current prime
   long long int  low_value;    // Lowest value on this proc
   char           *marked;      // Portion of 2,...,'n'
   long long int  n;            // Sieving from 2, ..., 'n'
   int            p;            // Number of processes
   long long int  proc0_size;   // Size of proc 0's subarray
   long long int  prime;        // Current prime
   long long int  size;         // Elements in 'marked'

   MPI_Init (&argc, &argv);

   // Start the timer

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   // Size of Process sero. This must have all prime roots that
   // will multiply through n
   // therefore it must be smaller than sqrt(n)
   // This ratio is critical to the indexing below.
   // and if changing the index, then change this ratio.
   proc0_size = (n-1)/p;
   // check to make sure all prime roots are within process zero
   if ((2 + proc0_size) < (long long int) sqrt((double) n)) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      exit (1);
   }

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

   // These are conditions that will be used to guide the algorithm
   // They will not neccessarilly be used in the algorithm calculations.
   low_value = 2 + id * proc0_size; //does not need to be integer
   high_value = 1 + (id+1) * proc0_size; //same, does not need to be int.
   size = high_value - low_value + 1; // This will be an integer.

   /* Allocate this process's share of the array. */
   marked = (char *) malloc (size);

   if (marked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit (1);
   }
   
   // initializing the algorithm.
   for (i = 0; i < size; i++) marked[i] = 0; // marks all zero to start
   // id is equal to 0
   if (!id){ 
      index = 0; // set index at beginning to start the algorithm
//      printf("when !id before do loop\n");
//      printf ("id is equal to: %d \n", id);
//      printf ("index is equal to: %d \n", index);
//      printf("-----------------------------------\n");
   }
   prime = 2;

//   printf("Right before do loop\n");
//   printf ("id is equal to: %d \n", id);
//   printf ("index is equal to: %d \n", index);
//   printf ("size is equal to: %d \n", size);
//   printf ("prime is equal to: %d \n", prime);
//   printf ("low_value is equal to: %d \n", low_value);
//   printf ("high_value is equal to: %d \n", high_value);
//   printf ("proc0_size is equal to: %d \n", proc0_size);
//   printf ("first is equal to: %lld \n", first);
//   printf("-----------------------------------\n\n\n");

   long long int times_through_loop = 0;

   // Start a while loop
   do {
      times_through_loop += 1;
//      printf("This is the %lld time through the loop\n\n\n", times_through_loop);
      // first is the base value to eliminate non-primes.
      // first is the first index of non-prime values, which has not been eliminated.
      if (prime * prime > low_value)
         first = prime * prime - low_value;
      else { // IF the remainder is zero (if devisable by...)
         if (!(low_value % prime)) first = 0;
         else first = prime - (low_value % prime);
      }
      // i is the local index
      // this loop eliminates all multiples of "prime"
      for (i = first; i < size; i += prime){ 
         marked[i] = 1;
//         printf("i is equal to %lld \n", i);
//         printf("the indexed value of i should be %lld \n", i+low_value);
      }
//      printf("\n\n\n");
      if (!id) {
         while (marked[++index]);
         prime = index + 2;
      }
//      printf("INSIDE DO LOOP\n");
//      printf ("id is equal to: %d \n", id);
//      printf ("index is equal to: %d \n", index);
//      printf ("size is equal to: %d \n", size);
//      printf ("prime is equal to: %d \n", prime);
//      printf ("low_value is equal to: %d \n", low_value);
//      printf ("high_value is equal to: %d \n", high_value);
//      printf ("proc0_size is equal to: %d \n", proc0_size);
//      printf ("first is equal to: %lld \n", first);
//      printf("-----------------------------------\n\n");
      if (p > 1) MPI_Bcast (&prime,  1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
   } while (prime * prime <= n);
   // END while loop

   count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;
   if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_LONG_LONG_INT, MPI_SUM,
      0, MPI_COMM_WORLD);

   /* Stop the timer */

   elapsed_time += MPI_Wtime();


   /* Print the results */

   if (!id) {
      printf ("There are %lld primes less than or equal to %lld\n",
         global_count, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   MPI_Finalize ();
   return 0;
}


