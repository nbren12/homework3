SRC = omp_solved2.c omp_solved3.c omp_solved4.c omp_solved5.c omp_solved6.c
EXE = $(SRC:.c=)


# CC     = gcc-4.9
CFLAGS += -fopenmp -O3
CFLAGS += -lm

all: bugs serial parallel

bugs: $(EXE)

serial: jacobi gs

parallel: jacobi-omp gs-omp


clean:
	rm -f jacobi jacobi-omp gs gs-omp $(EXE)  
