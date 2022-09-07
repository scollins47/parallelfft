Performance analysis on parallel variations of the FFT, made
by Jacob Hataway, Sammy Collins, and Cooper Steinour.

Our testing was done primarily using the file: sin_wav_10min_10000hz.txt.
This is a fairly large file (407M), and currently we do not have permission
to copy this to ./scratch/shared. Due to this, in order to test our code
the file can be found at /scratch/parallelfft/sin_wav_10min_10000hz.txt

There is additional input data in ./data, however, to run our implementations
with large sample sizes there needs to be enough samples in the input file.
Using a file with a smaller amount of samples than specified in the command
line will cause a segfault. For this reason it is best to use sin_wav_10min_10000hz,
which is sufficient enough to run up to 8,388,608 samples.

All of our source code can be compiled with the Makefile. All
implementations require powers of 2 for the number of samples.
To run use the following command formats:

FHT:
    OMP_NUM_THREADS=[n] ./fht [number of samples] [filename]
    Example:
        OMP_NUM_THREADS=2 ./fht 8388608 ./sin_wav_10min_10000hz.txt

dit_rad2:
    OMP_NUM_THREADS=[n] ./dit_rad2 [number of samples] [filename]
    Example:
        OMP_NUM_THREADS=2 ./dit_rad2 8388608 ./sin_wav_10min_10000hz.txt

dif_rad2:
    OMP_NUM_THREADS=[n] ./dif_rad2 [filename] [number of samples] [radix]
    Note that [radix] should always be 2, as radix 4 was never parallelized.
    Example:
        OMP_NUM_THREADS=2 ./dif_rad2 ./sin_wav_10min_10000hz.txt 8388608 2