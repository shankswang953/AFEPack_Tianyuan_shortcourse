CXX ?= c++

AFEPACK_PREFIX ?= $(HOME)
OPENBLAS_PREFIX ?= /opt/local
BOOST_INCLUDE ?= /opt/local/libexec/boost/1.81/include

CPPFLAGS := -I$(AFEPACK_PREFIX)/include \
            -I$(OPENBLAS_PREFIX)/include/openblas \
            -I$(BOOST_INCLUDE)
CXXFLAGS := -O2 -g -std=c++20 -Wall -Wextra
LDFLAGS := -L$(AFEPACK_PREFIX)/lib \
           -L$(OPENBLAS_PREFIX)/lib \
           -Wl,-rpath,$(OPENBLAS_PREFIX)/lib
LDLIBS := -lAFEPack -lopenblas

.PHONY: all clean run

all: shape_optimization

shape_optimization: shape_optimization.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

shape_optimization.o: shape_optimization.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

run: shape_optimization
	./shape_optimization output 4

clean:
	rm -f shape_optimization shape_optimization.o
