#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "timer.h"
/***********************************************************
* static void fht_dit_rec(double *x, unsigned long n, int nbranch)
*
* Purpose:
*   Computes the discrete Hartley transform of a real
*   sequence x[0..n-1], using a radix 2 decimation-in-time
*   FHT.  If the computation is small enough to fit in
*   cache, it is done iteratively.  Otherwise, it is done
*   recursively until the recursion descends to cache-sized
*   chunks.
*
*   n must be a power of 2.  Entering this function, x[]
*   must be in bit-reversed order.  On return, x[] contains
*   the Hartley transform, returned to normal order.
*
*   To support OpenMP parallelism, nbranch keeps track of
*   the number of active transforms at a given recursion
*   level. On the first call to this function, nbranch
*   should be 1.  It is then doubled for each recursion.
*
* Arguments:
*   double *x       - array of n doubles, representing n
*                     real numbers
*   unsigned long n - dimension of x, must be a power of 2
*   int nbranch     - number of transforms at this recursion
*                     level
************************************************************/
const double PI = 3.1415926535897932;
static void fht_dit_rec(double *x, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long j, jmax, k, nh, nq;

    if (n == 1)
        return;
    nh = n >> 1;
    nq = nh >> 1;
    nbranch <<= 1;
    #if (_OPENMP >= 200203)
    if (nbranch > omp_get_max_threads()) {
        fht_dit_rec(x     , nh, nbranch);
        fht_dit_rec(x + nh, nh, nbranch);
    } else {
        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            fht_dit_rec(x     , nh, nbranch);
            #pragma omp section
            fht_dit_rec(x + nh, nh, nbranch);
        }
    }
    #else
    fht_dit_rec(x     , nh, nbranch);
    fht_dit_rec(x + nh, nh, nbranch);
    #endif
    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    jmax = nq + nh;
    c = 1.0;
    s = 0.0;
    for (j = nh + 1, k = n - 1; j < jmax; ++j, --k) {
        double tmp, u, v;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = u * c + v * s;
        x[k] = u * s - v * c;
    }

    for (j = 0; j < nh; ++j) {
        double u, v;
        u = x[j];
        v = x[j + nh];
        x[j] = u + v;
        x[j + nh] = u - v;
    }
    return;
}   /* fht_dit_rec() */

static long getSizeOfArray(char *fname){
    FILE *file = fopen(fname, "r");
    char line[256];
    memset(line, 0, 256);
    long numSamples;
    //first line is always the number of samples taken
    fgets(line, sizeof(line), file);
    numSamples = strtol(line, NULL, 10);
    fclose(file);
    return numSamples;
}

static double* toArray (char* fname){
    FILE* file = fopen(fname, "r");
    char line[256];
    memset(line, 0, 256);
    long numSamples;
    //first line is always the number of samples taken
    fgets(line, sizeof(line), file);
    numSamples = strtol(line, NULL, 10);
    double *sampleArray = (double *) malloc(numSamples * sizeof(double));
    
    memset(sampleArray, 0, numSamples);
    int index = 0;
    while (fgets(line, sizeof(line), file))
    {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        sampleArray[index++] = strtod(line, NULL);
        //printf("%s", line);
    }
    /* may check feof here to make a difference between eof and io failure -- network
       timeout for instance */

    fclose(file);
    return sampleArray;
}
void printArray(double*arr, int length, char* fname){
    if (fname){
        FILE* file = fopen(fname, "w");
        for (int i = 0; i < length; i++)
        {
            char temp[100];
            snprintf(temp, 99,"%f %d\n", arr[i], i);
            fputs(temp, file);
        }
    } else {
        for (int i = 0; i < length; i++)
        {
            printf("%f %d\n", arr[i], i);
        }
    }
}
// testing dimensions: 
// dimensions must be > sample size
/*
    sin_output.txt - 88,200
    africa_output.txt - 6,570,432 
    sinwav_10000hz - 28 million samples
    2^10 1024
    2^15 32768
    2^20 1048576
    2^22 4,194,304 - highest w/o seg fault for africa and sin

    2^25 33,554,432
    2^30 1073741824
    2^32 4294967296
*/
int main(int argc, char** argv){
    char *fname;
    unsigned long ndimensions;
    bool writeout;
    if (argc < 3){
        printf("USAGE: <number of dimensions> <filename> <writeout default -f alse>\n NUM DIMENSIONS MUST BE POWER OF 2 \n");
        exit(1);
    }
    ndimensions = atoi(argv[1]);
    fname = argv[2];
    if (argc == 4){
        writeout = true;
    } else {
        writeout = false;
    }
    long length = getSizeOfArray(fname);
    double *arr = toArray(fname);
    
    START_TIMER(timer)
    fht_dit_rec(arr, ndimensions, 2);
    STOP_TIMER(timer)
#ifdef _OPENMP
    printf("TIME: %8.4f seconds\n", GET_TIMER(timer));
    if (writeout)
        printArray(arr, length, "parallelOut.txt");
#else
    if (writeout)
        printArray(arr, length, "serialOut.txt");
    printf("TIME: %8.4f seconds\n", GET_TIMER(timer));
#endif

    free(arr);
}
