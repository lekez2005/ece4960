LINK.o = $(LINK.cc)
CFLAGS=-Wall
CXXFLAGS=-Wall

.DEFAULT_GOAL := all

main: main.o

main.o: base_matrix.h full_matrix.h full_matrix.cpp sparse_matrix.h sparse_matrix.cpp

.PHONY: all
all: main
	@true

.PHONY: clean
clean:
	rm -fr *.o
