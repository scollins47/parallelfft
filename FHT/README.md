# Performance Analysis of parallel Fast Hartley Transformation
Performance analysis of Iterative Radix-2 Decimation in Time Fast Hartley Transformation

### Description of files:
#### **par_sammys.c**
    This is the parallel version of the Fast Hartley Transformation, must be compile with -fopenmp. To specify the number of threads `OMP_NUM_THREADS=X` where X is the desired number of threads. Syntax for running the main method in this file is as follows:
     `./par_sammys <desired dimension> <file to read samples from> <Left blank if just for timing, anything else to write to the file parallelOut.txt>` Specifying the number of threads using OMP_NUM_THREADS is not required, but lack of doing so will result in an unparallelized performance.

#### **sammys.c**
    This is the serial version unchanged in the implementation of the FHT, however, a main was added along with helper methods to read and write to files. Usage is as follows:
    `./sammys <# dimensions> <inputfile> <left blank if just for timing, anything else will write to the file serialOut.txt>`

#### Input file format
    The preloaded files can be quite large, and are as follows:
    sin_output.txt - Sampled 3 second 440 Hz sine wave
    africa_output.txt - Sampled Africa - Toto, (I wish I could say I licensed it but :\)
    sin_wav_10min_10000hz.txt - Sampled 10 minute 10000 Hz sine wave
    Input files should be as follows: 
    The first line should be the number of samples in the file (used for instantiating an array without having to deal with dynamic memory allocation)
    All lines following must have at least 1 number, that number being the sample, for graphing purposes many of the input files will be in (y,x) format however the x coordinate will not be read by the program, and is strictly for graphing purposes.

#### Header Files:
    timer.h - used for parallel and serial timing in both sammys.c and par_sammys.c

#### Output file format
    As is for the input file, the output file is formatted in a (y, x) format where y is the transformed sample, and x is used for strictly graphing purposes.


## How to use
    I have included a make file for ease of use, included in the make file are the following commands: make ____
    default - compiles and links both versions of the FHT, serial, and parallel with openmp.
    parallel - compiles and links just the parallel version of the FHT, subsequently linking openmp.
    serial - compiles and links just the serial version of the FHT.
    test - compares the output of the serial version, and the parallel version using diff, if correct there should be no output after the result of diff.
    testExhaustive - tests all possible combinations of the preloaded input files, and dimensions, writing them to testResults.txt formatted like so:
    -------- Version ----- ------ Num Dimensions ----- --- Num Samples -----
    Time 1
    Time 2
    Time 3
    -------- Version ----- ------ Num Dimensions ----- --- Num Samples -----
    Running the combination 3 times before advancing to the next one.