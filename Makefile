CC = gcc
CFLAGS = -W -Wall
LDFLAGS = -llapack -lm
EXEC = PRR

all: $(EXEC)

run: $(EXEC)
	./PRR

$(EXEC):
	$(CC) PRR.c -o $@ $(LDFLAGS)