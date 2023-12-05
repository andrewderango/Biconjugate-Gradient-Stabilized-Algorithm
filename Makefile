CC = gcc
CFLAGS = -Wall -Wextra -g -lm -O2 -fprofile-arcs -ftest-coverage $(shell pkg-config --cflags libpng)
LIBS = $(shell pkg-config --libs libpng)

all: bicgstab

bicgstab: main.c functions.c functions.h
	$(CC) $(CFLAGS) -o bicgstab main.c functions.c $(LIBS)

clean:
	rm -f bicgstab
	rm -rf bicgstab.dSYM
	rm -f *.gcda
	rm -f *.gcno
	rm -f *.gcov