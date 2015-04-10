SRC = omp_solved2.c
EXE = $(SRC:.c=)


CC     = gcc
CFLAGS = -fopenmp

all: $(EXE)
