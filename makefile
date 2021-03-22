CC=gcc
CFLAGS=-I.

strassen: strassen.c pcg_basic.c pcg_basic.h
	$(CC) $(CFLAGS) strassen.c pcg_basic.c -o randmst -O3