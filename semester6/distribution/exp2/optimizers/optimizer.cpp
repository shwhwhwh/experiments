#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <iostream>

#define LL long long
#define ODD(a) (a * 2 - 1)
LL idx;

LL find_primes(LL n, char *mark_mat, int p, int id, bool step) {
    LL count = 0;        /* Local non-prime count */
    LL first;        /* Index of first multiple */
    LL high_value;   /* Highest value on this proc */
    LL low_value;    /* Lowest value on this proc */
    LL i;
    LL index = 0;        /* Index of current prime */
    LL prime;        /* Current prime */
    LL size;         /* Elements in 'marked' */
    LL n_div_2;        /* n divided by 2 */
    LL head;
    char *marked;       /* Portion of 2,...,'n' */

    //data preprocess
    n_div_2 = (n + 1) / 2;
    low_value = 2 + id * (n_div_2 - 1) / p;
    high_value = 1 + (id + 1) * (n_div_2 - 1) / p;
    size = high_value - low_value + 1;
    low_value = ODD(low_value);
    if (!id && step)low_value = ((LL) sqrt((double) n) + 2) / 2 * 2 + 1;
    high_value = ODD(high_value);
    head = ((low_value - 3) / 2) * (!step);
    if (!id && step)head = (low_value - 3) / 2;
    //finding of primes
    if (!step) {
        prime = 3;
        do {
            if (prime * prime > low_value)
                first = (prime * prime - low_value) / 2 + head;
            else if (!(low_value % prime))first = head;
            else if (!((prime - (low_value % prime)) % 2))first = (prime - (low_value % prime)) / 2 + head;
            else first = (prime * 2 - (low_value % prime)) / 2 + head;
            for (i = first; i < size + head; i += prime) {
                if (!mark_mat[i])++count;
                mark_mat[i] = 1;
            }
            if (!id) {
                while (mark_mat[++index]);
                prime = index * 2 + 3;
            }
            if (p > 1) MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
        } while (prime * prime <= n);
    } else {
        if (id) {
            /* Allocate this process's share of the array. */
            marked = (char *) malloc(size);
            if (marked == NULL) {
                printf("Cannot allocate enough memory\n");
                MPI_Finalize();
                exit(1);
            }
            memset(marked, 0, size);
        } else marked = mark_mat;
        prime = 3;
        do {
            if (prime * prime > low_value)
                first = (prime * prime - low_value) / 2 + head;
            else if (!(low_value % prime))first = head;
            else if (!((prime - (low_value % prime)) % 2))first = (prime - (low_value % prime)) / 2 + head;
            else first = (prime * 2 - (low_value % prime)) / 2 + head;
            for (i = first; i < size; i += prime) {
                if (!marked[i])++count;
                marked[i] = 1;
            }
            while (mark_mat[++index]);
            prime = index * 2 + 3;
        } while (prime * prime <= n);
    }
    return count;
}

int main(int argc, char *argv[]) {
    double elapsed_time; /* Parallel execution time */
    LL global_count = 0; /* Global non-prime count in odds from 5 to n*/
    LL count = 0;        /* Local non-prime count */
    int id;           /* Process ID number */
    char *marked_0;
    LL n;            /* Sieving from 2, ..., 'n' */
    int p;            /* Number of processes */
    LL proc0_size;   /* Size of proc 0's subarray */
    int shm_id;         /* id for shared memory of shared prime */
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
    const LL SIZE = (LL) sqrt((double) n) + 1;
    n_div_2 = (n + 1) / 2;

    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = (n - 1) / p;
    if ((2 + proc0_size) < (int) sqrt((double) n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate shared memory of SIZE which is used to store primes in 3 to sqrt(n) */
    if (!id)shm_id = shmget(IPC_PRIVATE, (n_div_2 - 1) / p * sizeof(char), IPC_CREAT | IPC_EXCL | 0666);
    MPI_Bcast(&shm_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (shm_id < 0) {
        std::cerr << "shmget error" << std::endl;
        exit(1);
    }
    marked_0 = (char *) shmat(shm_id, NULL, 0);



    /* find primes from 3 to sqrt(n) + 1 */
    count += find_primes(SIZE, marked_0, p, id, false);
    shmdt(marked_0);
    MPI_Barrier(MPI_COMM_WORLD);

    /* find primes <= n */
    marked_0 = (char *) shmat(shm_id, NULL, 0);
    count += find_primes(n, marked_0, p, id, true);
    shmdt(marked_0);
    
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
               0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();
    if(!id)shmctl(shm_id, IPC_RMID, NULL);
    
    /* Print the results */

    if (!id) {
        printf("There are %lld primes less than or equal to %lld\n",
               (n + 1) / 2 - global_count - 1 + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}