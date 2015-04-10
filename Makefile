SRC = omp_solved2.c omp_solved3.c
EXE = $(SRC:.c=)


CC     = gcc
CFLAGS = -fopenmp

all: $(EXE)
