# https://stackoverflow.com/questions/25077300/include-source-files-using-path-in-makefile

SOURCES  := $(wildcard *.cpp)
EXECS:= $(SOURCES:.cpp=)

PYTHON_DIR = /home/lekez2005/anaconda2/envs/caffe2

MATRIX_LIB = -I../../lib

CXXFLAGS = $(MATRIX_LIB) -Wall -g
LDFLAGS = -lm -lpython2.7 -g

CXXFLAGS +=  -I/usr/include/python2.7

define CONDA
CXXFLAGS += -I$(PYTHON_DIR)/include/python2.7
CXXFLAGS += -DWITHOUT_NUMPY
# CXXFLAGS += -I$(PYTHON_DIR)/lib/python2.7/site-packages/numpy/core/include
LDFLAGS +=  -L$(PYTHON_DIR)/lib -Wl,-rpath,$(PYTHON_DIR)/lib
endef

CXXFLAGS += -std=c++11
LDFLAGS += -std=c++11

# need to put sources before flags to prevent python symbol errors
$(EXECS): $(SOURCES) 
	g++ $@.cpp -o $@ $(CXXFLAGS) $(LDFLAGS)

.PHONY: all

all: $(EXECS)

.PHONY: clean

clean:
	rm $(EXECS)
