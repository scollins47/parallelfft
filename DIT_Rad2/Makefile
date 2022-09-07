CC = gcc
CFLAGS= -std=c99 -Wall -g -lm

default: old test2 test2_serial

old: old.c timer.h
		$(CC) $(CFLAGS) -o old old.c -fopenmp

test2: test2.c timer.h
		$(CC) $(CFLAGS) -o test2 test2.c -fopenmp

test2_serial: test2.c timer.h
		$(CC) $(CFLAGS) -O2 -Wno-unknown-pragmas -o test2_serial test2.c

clean:
	rm -f old test2 test2_serial

