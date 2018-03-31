# https://stackoverflow.com/questions/25077300/include-source-files-using-path-in-makefile

SOURCES  := $(wildcard *.cpp)
EXECS:= $(SOURCES:.cpp=)

MATRIX_LIB = -I../../lib

CXXFLAGS = $(MATRIX_LIB) -Wall -g
LDFLAGS = -lm -g

.PHONY: all

all: $(EXECS)

.PHONY: clean

clean:
	rm $(EXECS)
