#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "timer.h"

#define TRSIZ 8
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
/******************************************************************************
* FUNCTION: diftt()
*
* PURPOSE:
* - Calculates Fast Fourier Transform of given data series using bit
* reversal after the FFT. This algorithm follows the butterfly diagram
* of the radix-2 Decimation In Frequency (DIF) implementation of FFT,
* which is modified from the previous Decimation In Time (DIT) algorithm.
* The function replaces the (input) array
* data[1..2*nn] by its discrete Fourier transfor, if the parameter isign
* is given as -1; or replaces data[1..2*nn] by nn times its inverse
* discrete Fourier transform, if isign is given as 1. The array of data[]
* is a complex array of length nn or, equivalently, a real array of
* length 2*nn. nn MUST be an integer power of 2, but this is not checked!
*
* The needed increments of trigonometric functions are handled via
* a clever identity using the fact that:
*
* exp(ia + ib) = exp(ia)*exp(ib) = exp(ia) + Z*exp(ia),
*
* where a is the starting angle and b is the increment. Solving for Z
* gives:
*
* Z = exp(ib) - 1 = -(1 - cos(b)) + isin(b) = -2*sin^2(b/2) + isin(b).
*
* Then the increments can be calculated with respect to the zero-angle
* set by wr (omega-real) and wi (omega-imaginary) by relation:
*
* (wr + wi) + Z*(wr + wi) = (wr + wi) + [-2*sin^2(b/2) + isin(b)](wr + wi).
*
* By setting wpr (omega-phase-real) = -2*sin^2(b/2) and wpi = sin(b),
* one has
*
* (wr + wi) + (wpr + wpi)*(wr + wi) = exp(ia) + Z*exp(ia),
*
* where the real part is: wr + wpr*wr - wpi*wi
* and the imaginary part is: wi + wpr*wi + wpi*wr.
* The actual incremental parts here are:
* wpr*wr - wpi*wi for real part and wpr*wi + wpi*wr for imaginary part.
*
*
* INPUT:
* - data[] = Array containing the data to be transformed
* The transformed data is stored back to this array
* so that the real and imaginary parts are following
* each other --> size of array = 2*size of data
* - nn = size of data = size of array/2, has to be a power of 2
* - isign = if -1, calculates the FFT;
* if 1, calculates the IFFT i.e. the inverse.
*
* OUTPUT: - (void)
*/
void diftt()
{
    // Initialize variables
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    int istep;
    int N = TRSIZ;
    int i = 0;
    int j = 0;
    int n = N*2;
    int  k = 0;
    int m = 0;
    int isign = -1;
    int mmax = n/2;

    double data1[2*TRSIZ] = {0,0,1,0,4,0,9,0,2,0,3,0,4,0,5,0};
    double *data;
    data = &data1[0] - 1;

    // calculate the FFT
    while (mmax >= 2) {
        istep = mmax << 1;
        theta = isign*(6.28318530717959/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = data[i];
                tempi = data[i+1];
                data[i] = data[i] + data[j];
                data[i+1] = data[i+1] + data[j+1];
                data[j] = tempr - data[j];
                data[j+1] = tempi - data[j+1];
                tempr = wr*data[j] - wi*data[j+1];
                tempi = wr*data[j+1] + wi*data[j];
                data[j] = tempr;
                data[j+1] = tempi;
                printf("\ni = %d ,j = %d, m = %d, wr = %f , wi = %f",(i-1)/2,(j-1)/2,m,wr,wi);
            }
            printf("\nm = %d ,istep = %d, mmax = %d, wr = %f , wi = %f, Z = %f"
            ,m,istep,mmax,wr,wi,atan(wi/wr)/(6.28318530717959/(1.0*n/2)));
            wtemp = wr;
            wr += wtemp*wpr - wi*wpi;
            wi += wtemp*wpi + wi*wpr;
        }
        mmax = mmax/2;
    }

    //DONT KNOW IF THIS IS NECESSARY YET????////
    // do the bit-reversal
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j +=m;
    }
    // print the results
    printf("\nFourier components from the DIF algorithm:");
    for (k = 0; k < 2*N; k +=2 )
        printf("\n%f %f", data[k+1], data[k+2]);
} // end of diftt()

// Main method to run the class
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






































// #define TRSIZ 8
// #define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
// /******************************************************************************
// * FUNCTION: dittt()
// *
// * PURPOSE:
// * - Calculates Fast Fourier Transform of given data series using bit
// * reversal prior to FFT. This function is directly taken from the book
// * "Numerical Recipes in C". The algorithm follows the butterfly diagram
// of the radix-2 Decimation In Time (DIT) implementation of FFT, where
// * the bit reversal is mandatory. The function replaces the (input) array
// * data[1..2*nn] by its discrete Fourier transfor, if the parameter isign
// * is given as -1; or replaces data[1..2*nn] by nn times its inverse
// * discrete Fourier transform, if isign is given as 1. The array of data[]
// * is a complex array of length nn or, equivalently, a real array of
// * length 2*nn. nn MUST be an integer power of 2, but this is not checked!
// *
// * The needed increments of trigonometric functions are handled via
// * a clever identity using the fact that:
// *
// * exp(ia + ib) = exp(ia)*exp(ib) = exp(ia) + Z*exp(ia),
// *
// * where a is the starting angle and b is the increment. Solving for Z
// * gives:
// *
// * Z = exp(ib) - 1 = -(1 - cos(b)) + isin(b) = -2*sin^2(b/2) + isin(b).
// *
// * Then the increments can be calculated with respect to the zero-angle
// * set by wr (omega-real) and wi (omega-imaginary) by relation:
// *
// * (wr + wi) + Z*(wr + wi) = (wr + wi) + [-2*sin^2(b/2) + isin(b)](wr + wi).
// *
// * By setting wpr (omega-phase-real) = -2*sin^2(b/2) and wpi = sin(b),
// * one has
// *
// * (wr + wi) + (wpr + wpi)*(wr + wi) = exp(ia) + Z*exp(ia),
// *
// * where the real part is: wr + wpr*wr - wpi*wi
// * and the imaginary part is: wi + wpr*wi + wpi*wr.
// * The actual incremental parts here are:
// * wpr*wr - wpi*wi for real part and wpr*wi + wpi*wr for imaginary part.
// *
// *
// * INPUT:
// * - data[] = Array containing the data to be transformed
// * The transformed data is stored back to this array
// * so that the real and imaginary parts are following
// * each other --> size of array = 2*size of data
// * - nn = size of data = size of array/2, has to be a power of 2
// * - isign = if -1, calculates the FFT;
// * if 1, calculates the IFFT i.e. the inverse.
// *
// * OUTPUT: - (void)
// */
// void dittt()
// {
//     // Initialize variables
//     double wtemp, wr, wpr, wpi, wi, theta;
//     double tempr, tempi;
//     int istep, mmax;
//     int N = TRSIZ;
//     int i = 0;
//     int j = 1;
//     int n = N*2;
//     int k = 0;
//     int m = 0;
//     int isign = -1;
//     double data1[2*TRSIZ] = {0,0,1,0,4,0,9,0,2,0,3,0,4,0,5,0};
//     double *data;
//     data = &data1[0] - 1;

//     // !!!!!May or may not need this !!!!!!!///
//     // do the bit-reversal
//     for (i = 1; i < n; i += 2) {
//     if (j > i) {
//     SWAP(data[j], data[i]);
//     SWAP(data[j+1], data[i+1]);
//     }
//     m = n >> 1;
//     while (m >= 2 && j > m) {
//     j -= m;
//     m >>= 1;
//     }
//     j +=m;
//     }

//     // calculate the FFT
//     mmax = 2;
//     while (n > mmax) {
//         istep = mmax << 1;
//         theta = isign*(6.28318530717959/mmax);
//         wtemp = sin(0.5*theta);
//         wpr = -2.0*wtemp*wtemp;
//         wpi = sin(theta);
//         wr = 1.0;
//         wi = 0.0;
//         for (m = 1; m < mmax; m += 2) {
//             for (i = m; i <= n; i += istep) {
//                 j = i + mmax;
//                 tempr = wr*data[j] - wi*data[j+1];
//                 tempi = wr*data[j+1] + wi*data[j];
//                 data[j] = data[i] - tempr;
//                 data[j+1] = data[i+1] - tempi;
//                 data[i] = data[i] + tempr;
//                 data[i+1] = data[i+1] + tempi;
//                 printf("\ni = %d ,j = %d, m = %d, wr = %f , wi = %f",(i-1)/2,(j-1)/2,m,wr,wi);
//             }
//             printf("\nm = %d ,istep = %d, mmax = %d, wr = %f , wi = %f, Z = %f"
//             ,m,istep,mmax,wr,wi,atan(wi/wr)/(6.28318530717959/(1.0*n/2)));
//             wtemp = wr;
//             wr += wtemp*wpr - wi*wpi;
//             wi += wtemp*wpi + wi*wpr;
//         }
//         mmax = istep;
//     }
//     // print the results
//     printf("\nFourier components from the DIT algorithm:");
//     for (k = 0; k < 2*N; k +=2 )
//     printf("\n%f %f", data[k+1], data[k+2]);
// } // end of dittt()

