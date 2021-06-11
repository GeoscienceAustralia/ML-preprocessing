#!/bin/env sh

CC=cc
all: core_c
clean:
	rm -f core_c.so	
core_c: 
	${CC} -fno-builtin -fPIC -shared -O3 -o core_c.so core_c.c -fopenmp -lgomp -Wall
