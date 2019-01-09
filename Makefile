CC = gcc
CFLAGS = -W -Wall
LDFLAGS = -llapack
EXEC = PRR

all: $(EXEC)

run: $(EXEC)
	./PRR

$(EXEC):
	$(CC) PRR.c -o $@ $(LDFLAGS)