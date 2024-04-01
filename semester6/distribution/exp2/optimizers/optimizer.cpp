#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include <string.h>

#define MIN(a, b)  ((a)<(b)?(a):(b))
#define LL long long
#define ODD(a) (a * 2 - 1)

LL get_first(LL p, LL low_v, LL high_v) {
    LL f;
    if (p * p > low_v)
        f = (p * p - low_v) / 2;
    else if (!(low_v % p))f = 0;
    else if (!((p - (low_v % p)) % 2))f = (p - (low_v % p)) / 2;
    else f = (p * 2 - (low_v % p)) / 2;
    return f;
}

int main(int argc, char *argv[]) {
    LL count = 0;        /* Local non-prime count */
    double elapsed_time; /* Parallel execution time */
    LL first;        /* Index of first multiple */
    LL global_count = 0; /* Global non-prime count in odds from 5 to n*/
    LL high_value;   /* Highest value on this proc */
    LL low_value;    /* Lowest value on this proc */
    LL i;
    int id;           /* Process ID number */
    LL index = 0;        /* Index of current prime */
    char *marked;       /* Portion of odds in 3,...,'n' */
    char *shared_marked;       /* Portion of odds in 3,...,'n' in process 0*/
    LL n;            /* Sieving from 2, ..., 'n' */
    int p;            /* Number of processes */
    LL proc0_size;   /* Size of proc 0's subarray */
    LL prime;        /* Current prime */
    LL size;         /* Elements in 'marked' */
    LL size_0;
    LL n_div_2;        /* n divided by 2 */
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
    size_0 = (n_div_2 - 1) / p;
    low_value = ODD(low_value);
    high_value = ODD(high_value);

    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = (n - 1) / p;
    if ((2 + proc0_size) < (int) sqrt((double) n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */

    marked = (char *) malloc(size);
    shared_marked = (char *) malloc((n_div_2 - 1) / p);
    if (marked == NULL || shared_marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    memset(marked, 0, size);
    memset(shared_marked, 0, size_0);
    prime = 3;
    LL low_value_0 = 2;
    LL high_value_0 = 1 + size_0;
    LL first_0;
    do {
        first = get_first(prime, low_value, high_value);
        first_0 = get_first(prime, low_value_0, high_value_0);
        for (i = first; i < size; i += prime) {
            if (!marked[i])++count;
            marked[i] = 1;
        }
        for (i = first_0; i < size_0; i += prime) {
            shared_marked[i] = 1;
        }
        while (marked[++index]);
        prime = index * 2 + 3;
    } while (prime * prime <= n);
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
               0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();

    /* Print the results */

    if (!id) {
        printf("There are %lld primes less than or equal to %lld\n",
               (n + 1) / 2 - global_count - 1 + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}
