// Mercredi 07 mai 2025
// Trying to emulate in C + OpenMP the test_parallel_sleep.sage
// Timings of this C code remain constant contrarily to the original
// SageMath.

// gcc -fopenmp -o test_openmp_sleep test_openmp_sleep.c
// (bash shell:)
// time ./test_openmp_sleep for using max number of threads (cores?)
// export NCPUS=<whatever> to set specific number of threads

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // For usleep
#include <omp.h>

void foo(int a, double T) {
    usleep((unsigned int)(T * 1000000)); // Sleep for T seconds (converted to microseconds)
}

void bar(int numcalls, double T, int ncpus) {
    for (int i = 0; i < numcalls; i++) {
        #pragma omp parallel for num_threads(ncpus)
        for (int a = 0; a < ncpus; a++) {
            foo(a, T);
        }
    }
}

int main() {
    int ncpus;
    // Check if NCPUS is set
    char *env_ncpus = getenv("NCPUS");
    if (env_ncpus != NULL) {
        ncpus = atoi(env_ncpus); // Use NCPUS value
    } else {
        ncpus = omp_get_max_threads(); // Default to max available threads (number of cores)
    }

    printf("I am set up to use %d threads.\n", ncpus);

    // Set the number of threads for OpenMP
    omp_set_num_threads(ncpus);

    int numcalls = 100; // Number of iterations
    double T = 0.005; // Individual sleep time in seconds

    // Call bar function
    bar(numcalls, T, ncpus);

    return 0;
}


// keeping this just in case
// Print the number of cores used in a parallel region
/* #pragma omp parallel */
/*   { */
/* #pragma omp single */
/*     printf("Actually using %d threads (cores?).\n", omp_get_num_threads()); */
/*   } */

