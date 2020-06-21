LIBS := `root-config --libs`
INCS := `root-config --cflags`
CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.exe : main.o random.o libreria.o point3D.o
	$(CC) random.o libreria.o point3D.o main.o -o main.exe  ${INCS} ${LIBS}

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS) ${INCS} ${LIBS}

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

point3D.o : point3D.cpp point3D.h
	$(CC) -c point3D.cpp -o point3D.o $(CFLAGS)

libreria.o : libreria.cpp libreria.h
	$(CC) -c libreria.cpp -o libreria.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
