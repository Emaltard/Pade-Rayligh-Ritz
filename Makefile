CC = mpicc
CFLAGS = -W -Wall
LDFLAGS = -fopenmp -llapacke -lm -g
EXEC = PRR

all: $(EXEC)

PRR: PRR.c mmio.c
	@$(CC) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o
	@rm -f PRR

run: 
	mpirun -n 2 ./$(EXEC) mm bcsstk01.mtx

test: 
	mpirun -n 1 ./$(EXEC) txt test.txt
