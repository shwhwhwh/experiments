#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include <string.h> 
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define LL long long
#define ODD(a) (a * 2 + 1)

int main(int argc, char* argv[])
{
	LL    count=0;        /* Local non-prime count */
    double elapsed_time; /* Parallel execution time */
    LL    first;        /* Index of first multiple */
    LL    global_count=0; /* Global non-prime count in odds from 5 to n*/
    LL    high_value;   /* Highest value on this proc */
    LL    low_value;    /* Lowest value on this proc */
    LL    i;
    int    id;           /* Process ID number */
    LL    index;        /* Index of current prime */
	char* marked;       /* Portion of 2,...,'n' */
	LL    n;            /* Sieving from 2, ..., 'n' */
	int    p;            /* Number of processes */
    LL    proc0_size;   /* Size of proc 0's subarray */
    LL    prime;        /* Current prime */
    LL    size;         /* Elements in 'marked' */
	LL    n_div_2;		/* n divided by 2 */
    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]);
    n_div_2 = (n + 1) / 2;
    
     /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */
    
    
    low_value = 2 + id * (n_div_2 - 1) / p;
    high_value = 1 + (id + 1) * (n_div_2 - 1) / p;
    size = high_value - low_value + 1;
    low_value = ODD(low_value);
    high_value = ODD(high_value);
    
    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = (n - 1) / p;

    if ((2 + proc0_size) < (int)sqrt((double)n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

     /* Allocate this process's share of the array. */
    
    marked = (char*)malloc(size);
    
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    memset(marked, 0, size);
    if(!id)index = 0;
    prime = 3;
    do {
        if (prime * prime > low_value)
            first = (prime * prime - low_value) / 2;
        else if(!(low_value % prime))first = 0;
        else if(!((prime-(low_value % prime))%2))first = (prime-(low_value % prime))/2;
        else first = (prime * 2 - (low_value % prime))/2;
        for(i = first; i < size; i += prime){
            marked[i] = 1;
            count++;
        }
        if(!id){
            while(marked[++index]);
            prime = index * 2 + 5;
        }
        if (p > 1) MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }while(prime * prime <= n);
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
        0, MPI_COMM_WORLD);
    
    /* Stop the timer */

    elapsed_time += MPI_Wtime();
    
     /* Print the results */

    if (!id) {
        printf("There are %lld primes less than or equal to %lld\n",
            (n+1)/2 - global_count - 1 + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}