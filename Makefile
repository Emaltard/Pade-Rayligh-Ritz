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
	mpirun -n 1 ./$(EXEC) txt test.txt