LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
 
%.o : %.cpp
	g++ -Wall -c $< ${INCS}

main.exe: qmc1d.o
	g++ -O3 -Wall -o main.exe $^ ${LIBS}

clean:
	rm *.o main.exe
