CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.exe : main.o random.o libreria.o city.o path.o
	$(CC) random.o libreria.o city.o path.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)

path.o : path.cpp path.h
	$(CC) -c path.cpp -o path.o $(CFLAGS)

libreria.o : libreria.cpp libreria.h
	$(CC) -c libreria.cpp -o libreria.o $(CFLAGS)

clean :
	rm *.o main.exe
