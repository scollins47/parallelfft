// from : https://lweb.cfa.harvard.edu/~spaine/am/download/src/transform.c

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

#ifndef FHT_UNIT_STRIDE
    #define FHT_UNIT_STRIDE 0
#endif
#define L1_CACHE_BYTES 0x8000
// 16 digits right of the decimal
const double PI = 3.1415926535897932;

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

/***********************************************************
* static void fht_dif_iter(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real
*   sequence x[0..n-1], using an iterative radix 2
*   decimation-in-frequency FHT.  n must be a power of 2.
*   Entering this function, x[] is in normal order.  On
*   return, x[] contains the Hartley transform, stored in
*   bit-reversed order.
*
* Arguments:
*   double *x       - array of n doubles, representing n
*                     real numbers
*   unsigned long n - dimension of x, must be a power of 2
************************************************************/

static void fht_dif_iter(double *x, unsigned long n)
{
    unsigned long m;

    for (m = n; m > 1; m >>= 1) {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double *xp;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
        }
        c = 1.0;
        s = 0.0;
        for (j = 1, k = mh - 1; j < mq; ++j, --k) {
            double tmp;
            double *xj, *xk;
            xj = x + j + mh;
            xk = x + k + mh;
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            for (i = 0; i < n; i += m) {
                double u, v;
                u = xj[i];
                v = xk[i];
                xj[i] = u * c + v * s;
                xk[i] = u * s - v * c;
            }
        }
    }
    return;
}   /* fht_dif_iter() */

/***********************************************************
* static void fht_dif_iter_seq(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real
*   sequence x[0..n-1], using an iterative radix 2
*   decimation-in-frequency FHT.  n must be a power of 2.
*   Entering this function, x[] is in normal order.  On
*   return, x[] contains the Hartley transform, stored in
*   bit-reversed order.
*
*   The two inner loops of the FHT computation are ordered
*   to favor sequential memory access at the expense of
*   redundant trig computations.  See J. Arndt, "Algorithms
*   for Programmers," online at http://www.jjj.de/fxt/.
*
* Arguments:
*   double *x       - array of n doubles, representing n
*                     real numbers
*   unsigned long n - dimension of x, must be a power of 2
************************************************************/

static void fht_dif_iter_seq(double *x, unsigned long n)
{
    unsigned long m;

    for (m = n; m > 1; m >>= 1) {
        double a, b, t;
        unsigned long i, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double c, s;
            double *xp;
            unsigned long j, k;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
            xp += mh;
            c = 1.0;
            s = 0.0;
            for (j = 1, k = mh - 1; j < mq; ++j, --k) {
                double u, v, tmp;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
                u = xp[j];
                v = xp[k];
                xp[j] = u * c + v * s;
                xp[k] = u * s - v * c;
            }
        }
    }
    return;
}   /* fht_dif_iter_seq() */

/***********************************************************
* static void fht_dif_rec(double *x, unsigned long n, int nbranch)
*
* Purpose:
*   Computes the discrete Hartley transform of a real
*   sequence x[0..n-1], using a radix 2 decimation-in-
*   frequency FHT.  If the computation is small enough to
*   fit in cache, it is done iteratively.  Otherwise, it is
*   done recursively until the recursion descends to
*   cache-sized chunks.
*
*   n must be a power of 2.  Entering this function, x[]
*   must be in normal order.  On return, x[] contains the
*   Hartley transform, in bit-reversed order.
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

static void fht_dif_rec(double *x, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long j, jmax, k, nh, nq;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / sizeof(double))) {
        if (FHT_UNIT_STRIDE)
            fht_dif_iter_seq(x, n);
        else
            fht_dif_iter(x, n);
        return;
    }
    nh = n >> 1;
    nq = nh >> 1;
    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    // can be parallelized
    for (j = 0; j < nh; ++j) {
        double u, v;
        u = x[j];
        v = x[nh + j];
        x[j] = u + v;
        x[nh + j] = u - v;
    }
    c = 1.0;
    s = 0.0;
    jmax = nq + nh;
    for (j = nh + 1, k = n - 1; j < jmax; ++j, --k) {
        double u, v, tmp;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = u * c + v * s;
        x[k] = u * s - v * c;
    }
    nbranch <<= 1;
    #ifdef _OPENMP
    if (nbranch > omp_get_max_threads()) {
        fht_dif_rec(x     , nh, nbranch);
        fht_dif_rec(x + nh, nh, nbranch);
    } else {
        //pragma omp sections num_threads(2)
        {// #pragma omp section
            #pragma omp task
            fht_dif_rec(x     , nh, nbranch);
            //#pragma omp section
            #pragma omp task 
            fht_dif_rec(x + nh, nh, nbranch);
            #pragma omp taskwait
        }
    }
    #else
    fht_dif_rec(x     , nh, nbranch);
    fht_dif_rec(x + nh, nh, nbranch);
    #endif
    return;
}   /* fht_dif_rec() */

/*
    sin_output.txt - 88,200
    africa_output.txt - 6,570,432 
    sinwav_10000hz - 28 million samples
    2^10 1024
    2^15 32768
    2^20 1048576
    2^22 4,194,304
    2^23 8,388,608
    2^25 33,554,432
    2^30 1073741824
    2^32 4294967296
*/
int main(int argc, char** argv){
    char *fname;
    unsigned long ndimensions;
    bool writeout;
    if (argc < 3){
        printf("USAGE: <number of dimensions> <filename> <writeout default false>\n NUM DIMENSIONS MUST BE POWER OF 2 \n");
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
#pragma omp parallel shared(arr, ndimensions)
#pragma omp single nowait
    fht_dif_rec(arr, ndimensions, 1);
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
