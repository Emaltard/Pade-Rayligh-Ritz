CC = mpicc
CFLAGS = -W -Wall
LDFLAGS = -fopenmp -llapacke -lm 
EXEC = PRR

all: $(EXEC)

PRR: PRR.c mmio.c
	@$(CC) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o
	@rm -f PRR

run4: 
	mpirun -n 4 ./$(EXEC) mm bcsstm01.mtx

run2: 
	mpirun -n 2 ./$(EXEC) mm bcsstm01.mtx

run1: 
	mpirun -n 1 ./$(EXEC) mm bcsstm01.mtx