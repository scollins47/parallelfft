#!/bin/bash

# build
make

echo "

== Strong Scaling =="
for t in 1 2 4 8 16 32 64; do
    echo "


== Parallel DIT for $t of input size 8388608 =="
    OMP_NUM_THREADS=$t ./test2 8388608 ../sin_wav_10min_10000hz.txt

done

echo "

== Weak Scaling =="
    echo "


== Parallel DIT for 1 of input size 1,024=="
    OMP_NUM_THREADS=1 ./test2 1024 ../sin_wav_10min_10000hz.txt

    echo "


== Parallel DIT for 2 of input size 32,768=="
    OMP_NUM_THREADS=2 ./test2 32768 ../sin_wav_10min_10000hz.txt

    echo "


== Parallel DIT for 4 of input size 1,048,576=="
    OMP_NUM_THREADS=4 ./test2 1048576 ../sin_wav_10min_10000hz.txt

    echo "


== Parallel DIT for 8 of input size 4,194,304=="
    OMP_NUM_THREADS=8 ./test2 4194304 ../sin_wav_10min_10000hz.txt


echo "

== Serial DIT for 1 =="
    OMP_NUM_THREADS=1 ./test2_serial 4194304 ../sin_wav_10min_10000hz.txt



