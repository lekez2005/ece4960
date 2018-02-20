LINK.o = $(LINK.cc)
CFLAGS=-Wall
CXXFLAGS=-Wall

main: main.o

main.o: base_matrix.h base_matrix.cpp full_matrix.h full_matrix.cpp sparse_matrix.h sparse_matrix.cpp

all: main

.PHONY: clean
clean:
	rm -fr *.o
