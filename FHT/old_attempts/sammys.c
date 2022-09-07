// fht_dif_iter() from : https://lweb.cfa.harvard.edu/~spaine/am/download/src/transform.c

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
// 16 digits right of the decimal
const double PI = 3.1415926535897932;
/***********************************************************
* static void fht_dit_iter(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real
*   sequence x[0..n-1], using an iterative radix 2
*   decimation-in-time FHT.  n must be a power of 2.
*   Entering this function, x[] must be in bit-reversed
*   order.  On return, x[] contains the Hartley transform,
*   returned to normal order.
*
* Arguments:
*   double *x       - array of n doubles, representing n
*                     real numbers
*   unsigned long n - dimension of x, must be a power of 2
************************************************************/

static void fht_dit_iter(double *x, unsigned long n)
{
    unsigned long m;
    for (m = 2; m <= n; m <<= 1) {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
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
    }
    return;
}   /* fht_dit_iter() */

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
/*
    sin_output.txt - 88,200
    africa_output.txt - 6,570,432 
    sinwav_10000hz - 28 million samples
    2^10 1024
    2^15 32768
    2^20 1048576
    2^22 4194304 - highest w/o seg fault

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
    fht_dit_iter(arr, ndimensions);
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
