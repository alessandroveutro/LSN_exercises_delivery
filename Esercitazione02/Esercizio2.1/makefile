CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.exe : main.o random.o libreria.o funzioni.h
	$(CC) random.o libreria.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

funzioni.o : libreria.cpp libreria.h
	$(CC) -c libreria.cpp -o libreria.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
