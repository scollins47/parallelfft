# Performance Analysis of parallel Decimation in Time Rad-2 FFT
Performance analysis of Iterative Radix-2 Decimation in Time FFT

### Description of files:
#### **test2.c**
    This is the parallel version of the implementation that can also be run serially, must be compile with -fopenmp. To specify the number of threads `OMP_NUM_THREADS=X` where X is the desired number of threads. Syntax for running the main method in this file is as follows:
     `OMP_NUM_THREADS=X ./test2 <desired dimension> <file to read samples from> 

#### Input file format
    The preloaded files can be quite large, and are as follows:
    sin_output.txt - Sampled 3 second 440 Hz sine wave
    africa_output.txt - Sampled Africa - Toto
    sin_wav_10min_10000hz.txt - Sampled 10 minute 10000 Hz sine wave
    The first line should be the number of samples in the file (used for instantiating an array without having to deal with dynamic memory allocation)
    All lines following must have at least 1 number, that number being the sample

#### Header Files:
    timer.h - used for parallel and serial timing in test2.c

#### Output file format
    The output file shows the transformed values for a serial program and then for the parallel program.

## How to use
    I have included a slurm file that will run both strong scaling and weak scaling on the sin_wav_10min_10000hz.txt which will run 
    by typing ./test_dif.sh


