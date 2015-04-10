SRC = omp_solved2.c omp_solved3.c omp_solved4.c omp_solved5.c
EXE = $(SRC:.c=)


CC     = gcc
CFLAGS = -fopenmp

all: $(EXE)
