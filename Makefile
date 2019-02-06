CC = mpicc
CFLAGS = -W -Wall
LDFLAGS = -llapacke -lm
EXEC = PRR

all: $(EXEC)

PRR: PRR.c mmio.c
	@$(CC) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o

run: 
	mpirun -n 1 ./$(EXEC) bcsstm01.mtx