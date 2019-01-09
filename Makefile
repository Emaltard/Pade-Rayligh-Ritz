CC=gcc
CFLAGS=-W -Wall
LDFLAGS= -llapack
EXEC=PRR

all: $(EXEC)

PRR: PRR.c
	@$(CC) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o