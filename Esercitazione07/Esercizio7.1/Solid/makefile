CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.exe : main.o libreria.o
	$(CC) libreria.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

libreria.o : libreria.cpp libreria.h
	$(CC) -c libreria.cpp -o libreria.o $(CFLAGS)

clean :
	rm *.o main.exe
