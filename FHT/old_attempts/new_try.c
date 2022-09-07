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
const double PI = 3.1415926535897932;
/**
    H is a bit reversed version of the signal
*/
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
// double reverseBits(double num){
//     double count = sizeof(num) * 8 - 1;
//     double reverse = num;
//     num >>= 1;
//     while (num) {
//         reverse <<= 1;
//         reverse |= num & 1;
//         num >>= 1;
//         count--;
//     }
//     reverse <<= count;
//     return reverse;
// }
// void reverseNums (double *arr, long length){
//     unsigned int numBits = sizeof(double) * 64;
//     for (long i = 0; i < length; i++){
//         arr[i] = reverseBits(arr[i]);
//     }
// }

void fht (double * H, long N){
    int P = 1; // where P is the number of threads
    int p = 0; // is p the current thread?
    long n = log(N) / log(2); // N points in the signal (num elements)
#ifdef _OPENMP
    P = omp_get_num_threads();
    p = omp_get_thread_num();
    printf(" RANK: %d\n", p);
#endif
    // for all (0 <= p <= P-1)
    long start = (p * N / (2 * P));
    long upto = ((p + 1) * N / (2 * P));
#pragma omp for
    for (long j = start; j < upto; j++) {
        double temp = H[2 * j + 1];
        H[2 * j + 1] = H[2 * j] - temp;
        H[2 * j] = H[2 * j] + temp;
    }

    // barrier synchronization
#pragma omp barrier
    for (long r = 1; r < n; r++ ){
        #pragma omp barrier
        // barrier synchronization
        if (r <= 2){
            // for all (0 <= p <= P-1)
            long start = p * N / (P * pow(2, r+1));
            long upto = (p + 1) * N/(P * pow(2, r+1));
#pragma     omp for
            for (int i = start; i < upto; i++){
                // T2 butterfly.
                long p2 = i*pow(2,r+1);
                long q2 = p2 + pow(2,r);
                long r2 = p2 + pow(2,r-1);
                long s2 = q2 + pow(2,r-1);
                double tmp1 = H[q2];
                double tmp2 = H[s2];
                H[q2] = H[p2] - tmp1;
                H[s2] = H[r2] - tmp2;
                H[p2] = H[p2] + tmp1;
                H[r2] = H[r2] + tmp2;

                for (long j = 1; j < pow(2, r-1); j++) {
                    // T1 butterfly
                    long p1 = p2 + j;
                    long q1 = p1 + pow(2,r);
                    long r1 = p2 + pow(2,r) - j;
                    long s1 = r1 + pow(2,r);
                    double Ci = cos((2 * PI * i) / N);
                    double Si = sin((2 * PI * i) / N);
                    double Cj = cos ((2 * PI * j) / N);
                    double Sj = sin ((2 * PI * j) / N);
                    double tmp1 = Ci*H[q1] + Si*H[s1];
                    double tmp2 = Cj*H[s1] + Sj*H[q1];
                    H[q1] = H[p1] - tmp1;
                    H[s1] = H[r1] - tmp2;
                    H[p1] = H[p1] + tmp1;
                    H[r1] = H[r1] + tmp2;
                }
            }
            
        } else {
            for (long i = 0; i < N / pow(2, r+1); i++){
                //compute t2 butterfly
                long p2 = i*pow(2,r+1);
                long q2 = p2 + pow(2,r);
                long r2 = p2 + pow(2,r-1);
                long s2 = q2 + pow(2,r-1);
                double tmp1 = H[q2];
                double tmp2 = H[s2];
                H[q2] = H[p2] - tmp1;
                H[s2] = H[r2] - tmp2;
                H[p2] = H[p2] + tmp1;
                H[r2] = H[r2] + tmp2;
                // all threads should calculate their own start and upto
                long start = p * pow(2, r - 1);
                long upto = (p + 1) * pow(2, r-1);
#pragma         omp for
                for (long j = start; j < upto; j++){
                    //compute t1 butterfly
                    long p1 = p2 + j;
                    long q1 = p1 + pow(2,r);
                    long r1 = p2 + pow(2,r) - j;
                    long s1 = r1 + pow(2,r);
                    double Ci = cos((2 * PI * i) / N);
                    double Si = sin((2 * PI * i) / N);
                    double Cj = cos ((2 * PI * j) / N);
                    double Sj = sin ((2 * PI * j) / N);
                    double tmp1 = Ci*H[q1] + Si*H[s1];
                    double tmp2 = Cj*H[s1] + Sj*H[q1]; // seg faults
                    H[q1] = H[p1] - tmp1;
                    H[s1] = H[r1] - tmp2;
                    H[p1] = H[p1] + tmp1;
                    H[r1] = H[r1] + tmp2;
                }
            }
        }
    }
}

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
    //reverseNums(arr, length);
    START_TIMER(timer)
#pragma omp parallel default(none) shared(arr, ndimensions)
    fht(arr, ndimensions);
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
