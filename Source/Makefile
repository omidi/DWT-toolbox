CC=g++
CFLAGS=-c -Wall -g3 -O3 -w
LDLIBS=-I boost_1_41_0/
LDFLAGS=
COMMON_SOURCES=Score.cpp Alignment.cpp BackgroundModel.cpp constants.cpp Decomposition.cpp Motif.cpp Window.cpp ParameterFile.cpp 
TARGET_SOURCES=DWT_model.cpp
TEST_SOURCES=Positional_Dependency_Posterior.cpp
COMMON_OBJECTS=$(COMMON_SOURCES:.cpp=.o)
TARGET_OBJECTS=$(TARGET_SOURCES:.cpp=.o)
TEST_OBJECTS = $(TEST_SOURCES:.cpp=.o)
EXECUTABLE=DWT_TFBS_prediction
TEST_EXECUTABLE=positional_dependency_posterior
$(info Unzipping the boost library. )
BOOST := $(shell tar -xf boost_1_41_0.tar.gz)

$(info Now compiling the programs. )
.PHONY: initialization all target tests

all: initialization target tests motevo clean

target: $(EXECUTABLE)

tests: $(TEST_EXECUTABLE)

$(EXECUTABLE): $(COMMON_OBJECTS) $(TARGET_OBJECTS)
	$(CC) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(TEST_EXECUTABLE): $(COMMON_OBJECTS) $(TEST_OBJECTS)
	$(CC) $(LDFLAGS) $(LDLIBS) $^ -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(LDLIBS) $< -o $@

motevo: motevo.c
	g++ -w -O2 motevo.c -o motevo

clean:
	rm -rf *.o
