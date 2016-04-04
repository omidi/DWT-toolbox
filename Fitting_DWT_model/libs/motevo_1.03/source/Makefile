all: cleanall ufe motevo clean

ufe: ufe1.o ufe2.o
	g++ -O2 runUFE.o evomodel.o -o runUFE

ufe1.o: evomodel.cpp
	g++ -O2 -c evomodel.cpp  

ufe2.o: runUFE.cpp
	g++ -O2 -c runUFE.cpp


motevo: motevo.c
	g++ -O2 motevo.c -o motevo


cleanall: 
	rm -rf *.o motevo runUFE

clean: 
	rm -rf *.o


