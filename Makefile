CC = mpicc
CFLAGS = -W -Wall
LDFLAGS = -fopenmp -llapacke -lm 
EXEC = PRR

all: $(EXEC)

PRR: PRR.c
	@$(CC) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o

run: 
	mpirun -n 1 ./$(EXEC)