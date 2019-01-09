CC = gcc
CFLAGS = -W -Wall
LDFLAGS = -llapacke -lm
EXEC = PRR

all: $(EXEC)

PRR: PRR.c
	@$(CC) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o

run: 
	./$(EXEC)