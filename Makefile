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

run: 
	mpirun -n 1 ./$(EXEC) mm plat1919.mtx

run4: 
	mpirun -n 4 ./$(EXEC) mm bcsstk01.mtx

test: 
	mpirun -n 1 ./$(EXEC) txt test.txt

test2: 
	mpirun -n 2 ./$(EXEC) txt test.txt