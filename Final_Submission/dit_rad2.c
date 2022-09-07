#include <string.h>
#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "timer.h"

#define PI 3.14159265358979323846
typedef unsigned long complex cplx;

// Reverse the bits
cplx* bit_reverse(cplx array[], int N){
    cplx temp;
    int n_bits = log2(N - 1) + 1;

    #pragma omp parallel for
    for (int i=0; i<N; i++){

        int r = 0;
        for (int j=0; j<n_bits; j++){
            int bit_j = (i >> j) & 1;
            r = r | (bit_j << (n_bits-j-1));
        }

        if (i < r) {
            temp = array[i];
            array[i] = array[r];
            array[r] = temp;
        }
    }

    return array;
}


cplx* W_fill(int N){
    cplx* W_array = (cplx*) malloc(N/2 * sizeof(cplx));
#   pragma omp parallel for
    for (int i=0; i<N/2; i++)
        W_array[i] = cexp(- 2 * I * i * PI / N);

    return W_array;
}


void fft_iterate(cplx array[], cplx W_array[], int N){
    array = bit_reverse(array, N);

    for (int layer=1; layer<log2(N - 1) + 1; layer++) {

        int group_size = 1 << layer;
        int group_num = N >> layer;
        int group_size_half = group_size >> 1;

        #pragma omp parallel for
        for (int group=0; group<group_num; group++){

            int index_0_i = group * group_size;

            for (int index_2_i=0; index_2_i<group_size_half*group_num; index_2_i += group_num) {

                int index_1_i = index_0_i + group_size_half;

                cplx f = array[index_0_i];;
                cplx g = array[index_1_i] * W_array[index_2_i];

                array[index_0_i] = f + g;
                array[index_1_i] = f - g;

                index_0_i ++;
            }
        }
    }
}

// Read in elements of a file
static cplx* toArray (char* fname){
    FILE* file = fopen(fname, "r");
    char line[256];
    memset(line, 0, 256);
    long numSamples;
    //first line is always the number of samples taken
    fgets(line, sizeof(line), file);
    numSamples = strtol(line, NULL, 10);
    cplx *sampleArray = (cplx *) malloc(numSamples * sizeof(cplx));

    memset(sampleArray, 0, numSamples);
    int index = 0;
    while (fgets(line, sizeof(line), file))
    {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        sampleArray[index++] = strtod(line, NULL);
    }
    /* may check feof here to make a difference between eof and io failure -- network
       timeout for instance */

    fclose(file);
    return sampleArray;
}

void printArray(cplx *arr, int length, char* fname){
        FILE* file = fopen(fname, "w");
        for (int i = 1; i < length; i++)
        {
            char temp[100];
            snprintf(temp, 99,"%.3Lf\n", creall(arr[i]));
            fputs(temp, file);
        }

}


int main(int argc, char** argv) {

    char *fname;
    int num_threads;
    unsigned long length;

#   ifdef _OPENMP
    num_threads = omp_get_max_threads();
#   else
    num_threads = 1;
#   endif

    if (argc < 3){
        printf("USAGE: <number of dimensions> <filename> NUM DIMENSIONS MUST BE POWER OF 2\n");
        exit(1);
    }

    length = atoi(argv[1]);
    fname = argv[2];

//    long length = getSizeOfArray(fname);
    cplx *arr = toArray(fname);

    START_TIMER(timer)
    cplx* W = W_fill(length);
    fft_iterate(arr, W, length);
    STOP_TIMER(timer)

    // print results
    printf("\nNthreads=%2d  Time Taken=%8.5fs\n",
            num_threads, GET_TIMER(timer));

#   ifdef _OPENMP
    //printArray(arr, length, "parallelOut.txt");
#   else
    //printArray(arr, length, "serialOut.txt");
#   endif

    free(arr);                 // Free arrays from memory after use
    free(W);
}


