CC=g++
CFLAGS=-c -Wall -g3 -O3 -w
LDLIBS=-I ../../lib/boost_1_41_0/
LDFLAGS=
SOURCES=Score.cpp Alignment.cpp BackgroundModel.cpp constants.cpp Decomposition.cpp Motif.cpp Window.cpp ParameterFile.cpp Positional_Dependency_Posterior.cpp
OBJECTS=$(SOURCES:.cpp=.o)
PROG=positional_dependency_posterior

all: $(PROG)

.cpp.o:
	$(CC) $(CFLAGS) $(LDLIBS) -o $@ $<

$(PROG): $(OBJECTS) 
	$(CC) -Wall -g3 -O3 -w $(OBJECTS) -o $@

.PHONY: clean
clean:
	-rm $(OBJECTS)

