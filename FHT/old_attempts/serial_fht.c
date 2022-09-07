#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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

// double reverseBits(double num){
//     double count = sizeof(num) * 8 - 1;
//     double reverse = num;
//     unsigned long *point = &num;
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

void fht(double *H, int N){
    int n = log(N) / log(2);
    for(int i = 0; i < N/2; i++){
        int temp = H[2 * i+1];
        H[2 * i+1] = H[2 * i] - temp;
        H[2 * i] = H[2 * i] + temp;
    }
    for(int r = 1; r < n; r++) {
        for(int i = 0; i < N/pow(2,r+1); i++) {
            int p2 = i*pow(2,r+1);
            int q2 = p2 + pow(2,r);
            int r2 = p2 + pow(2,r-1);
            int s2 = q2 + pow(2,r-1);
            int tmp1 = H[q2];
            int tmp2 = H[s2];
            H[q2] = H[p2] - tmp1;
            H[s2] = H[r2] - tmp2;
            H[p2] = H[p2] + tmp1;
            H[r2] = H[r2] + tmp2;
            for(int j = 1; j < pow(2,r-1); j++) {
                int p1 = p2 + j;
                int q1 = p1 + pow(2,r);
                int r1 = p2 + pow(2,r) - j;
                int s1 = r1 + pow(2,r);
                double Ci = cos((2 * PI * i) / N);
                double Si = sin((2 * PI * i) / N);
                double Cj = cos ((2 * PI * j) / N);
                double Sj = sin ((2 * PI * j) / N);
                int tmp1 = Ci*H[q1] + Si*H[s1];
                int tmp2 = Cj*H[s1] + Sj*H[q1];
                H[q1] = H[p1] - tmp1;
                H[s1] = H[r1] - tmp2;
                H[p1] = H[p1] + tmp1;
                H[r1] = H[r1] + tmp2;
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
    fht(arr, ndimensions);
    STOP_TIMER(timer)
    if (writeout)
        printArray(arr, length, "serialOut.txt");
    printf("TIME: %8.4f seconds\n", GET_TIMER(timer));


    free(arr);
}
