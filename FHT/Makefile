CC = gcc
CFLAGS= -std=c99 -Wall -W -g -lm

default: timer.h
		$(CC) $(CFLAGS) -o par_sammys par_sammys.c -fopenmp
		$(CC) $(CFLAGS) -o sammys sammys.c -lm
		$(CC) $(CFLAGS) -o par_sammys_rec par_sammys_rec.c -fopenmp

parallel: timer.h
		$(CC) $(CFLAGS) -o par_new_try new_try.c -fopenmp

serial: timer.h
		$(CC) $(CFLAGS) -o fht fht.c -lm

recursive: timer.h
		$(CC) $(CFLAGS) -o par_fht fht.c -fopenmp

wav:
		$(CC) $(CFLAGS) example.c wavfile.c -o example

clean:
		rm -f fht.o fht par_fht par_fht.o
test:
		OMP_NUM_THREADS=8 srun -n 1 par_fht 4194304 ../sin_wav_10min_10000hz.txt true
		srun -n 1 fht 4194304 ../sin_wav_10min_10000hz.txt true
		diff serialOut.txt parallelOut.txt
		OMP_NUM_THREADS=4 srun -n 1 par_fht 1024 ../sin_output.txt true
		srun -n 1 fht 1024 ../sin_output.txt true
		diff serialOut.txt parallelOut.txt

time:
	    srun -n 1 ./fht 8388608 ../sin_wav_10min_10000hz.txt
		OMP_NUM_THREADS=2 srun -n 1 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt
		OMP_NUM_THREADS=4 srun -n 1 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt
		OMP_NUM_THREADS=8 srun -n 1 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt
		OMP_NUM_THREADS=16 srun -n 1 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt
		OMP_NUM_THREADS=32 srun -n 1 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt

testExhaustive:
		echo "-----SERIAL----- ---- 1024 dimensions, 88,200 samples ----" >> testResults.txt
		./fht 1024 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 1024 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 1024 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----SERIAL----- ---- 32768 dimensions, 88,200 samples ----" >> testResults.txt
		./fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----SERIAL----- ---- 1024 dimensions, 6,570,432 samples ----" >> testResults.txt
		./fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----SERIAL----- ---- 32768 dimensions, 6,570,432 samples ----" >> testResults.txt
		./fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----SERIAL----- ---- 1048576 dimensions, 6,570,432 samples ----" >> testResults.txt
		./fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----SERIAL----- ---- 4194304 dimensions, 6,570,432 samples ----" >> testResults.txt
		./fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----SERIAL----- ---- 4194304 dimensions, 28,825,252 samples ----" >> testResults.txt
		./fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----SERIAL----- ---- 8388608 dimensions, 28,825,252 samples ----" >> testResults.txt
		./fht 8388608 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 8388608 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		./fht 8388608 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1024 dimensions, 88,200 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32768 dimensions, 88,200 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32768 dimensions, 88,200 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32768 dimensions, 88,200 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../sin_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1024 dimensions, 6,570,432 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1024 dimensions, 6,570,432 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1024 dimensions, 6,570,432 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1024 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32768 dimensions, 6,570,432 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32768 dimensions, 6,570,432 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32768 dimensions, 6,570,432 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1048576 dimensions, 6,570,432 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1048576 dimensions, 6,570,432 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1048576 dimensions, 6,570,432 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1048576 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 4194304 dimensions, 6,570,432 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 4194304 dimensions, 6,570,432 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 4194304 dimensions, 6,570,432 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 4194304 ../africa_output.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1024 dimensions, 28,825,252 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1024 dimensions, 28,825,252 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1024 dimensions, 28,825,252 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1024 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32,768 dimensions, 28,825,252 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32,768 dimensions, 28,825,252 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 32,768 dimensions, 28,825,252 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 32768 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1,048,576 dimensions, 28,825,252 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1,048,576 dimensions, 28,825,252 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 1,048,576 dimensions, 28,825,252 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 1048576 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 4,194,304 dimensions, 28,825,252 samples, 2 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=2 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 4,194,304 dimensions, 28,825,252 samples, 4 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=4 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 4,194,304 dimensions, 28,825,252 samples, 8 THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 4194304 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		echo "-----PARALLEL----- ---- 8,388,608 dimensions, 28,825,252 samples, 8 THREADS ----" >> testResults.txt
		THREADS ----" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt
		OMP_NUM_THREADS=8 ./par_fht 8388608 ../sin_wav_10min_10000hz.txt >> testResults.txt
		echo "" >> testResults.txt


