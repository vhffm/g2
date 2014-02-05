CC=gcc -Wall -ansi -pedantic

all: modpart_ssic

modpart_ssic: modpart_ssic.o ssio.o
	$(CC) -o modpart_ssic modpart_ssic.o ssio.o 

ssio.o: ssio.c
	$(CC) -c ssio.c

modpart_ssic.o: modpart_ssic.c
	$(CC) -c modpart_ssic.c

modpart_ssic.o: ssio.o ssio.h

clean:
	rm  *.o
