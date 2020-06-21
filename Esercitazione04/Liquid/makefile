LIBS := `root-config --libs`
INCS := `root-config --cflags`
CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.exe : MolDyn_NVE.o
	$(CC) MolDyn_NVE.o -o main.exe  ${INCS} ${LIBS}

MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)

clean :
	rm *.o main.exe

restart:
	cp start/config.0 .
