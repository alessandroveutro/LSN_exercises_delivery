CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.exe : Variational_MC.o random.o 
	$(CC) random.o Variational_MC.o -o main.exe

Variational_MC.o : Variational_MC.cpp Variational_MC.h
	$(CC) -c Variational_MC.cpp -o Variational_MC.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe
