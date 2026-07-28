CXX ?= c++

AFEPACK_PREFIX ?= $(HOME)
OPENBLAS_PREFIX ?= /opt/local
BOOST_INCLUDE ?= /opt/local/libexec/boost/1.81/include

CPPFLAGS := -I$(AFEPACK_PREFIX)/include \
            -I$(OPENBLAS_PREFIX)/include/openblas \
            -I$(BOOST_INCLUDE)
CXXFLAGS := -O2 -g -std=c++20
LDFLAGS := -L$(AFEPACK_PREFIX)/lib \
           -L$(OPENBLAS_PREFIX)/lib \
           -Wl,-rpath,$(OPENBLAS_PREFIX)/lib
LDLIBS := -lAFEPack -lopenblas

.PHONY: all clean

all: main

main: boundary_flux_dual.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

boundary_flux_dual.o: boundary_flux_dual.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f main boundary_flux_dual.o
