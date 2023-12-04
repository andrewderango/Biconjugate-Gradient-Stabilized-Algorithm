CC = gcc
CFLAGS = -Wall -Wextra -g -O2 $(shell pkg-config --cflags libpng)
LIBS = $(shell pkg-config --libs libpng)

all: bicgstab

bicgstab: main.o functions.o
	$(CC) $(CFLAGS) -o bicgstab main.o functions.o $(LIBS)

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c main.c

functions.o: functions.c functions.h
	$(CC) $(CFLAGS) -c functions.c

clean:
	rm -f *.o bicgstab
