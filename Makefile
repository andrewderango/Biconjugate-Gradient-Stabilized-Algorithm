CC = gcc
CFLAGS = -Wall -Wextra -g -lm -O2 $(shell pkg-config --cflags libpng)
LIBS = $(shell pkg-config --libs libpng)

all: bicgstab

bicgstab: main.c functions.c functions.h
	$(CC) $(CFLAGS) -o bicgstab main.c functions.c $(LIBS)

clean:
	rm -f bicgstab
	rm -rf bicgstab.dSYM