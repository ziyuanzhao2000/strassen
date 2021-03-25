CC=gcc
CFLAGS=-I.

strassen: strassen.c pcg_basic.c pcg_basic.h
	$(CC) $(CFLAGS) strassen.c -lm pcg_basic.c -o strassen -O3