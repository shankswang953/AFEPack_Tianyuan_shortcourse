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

.PHONY: all clean

all: generate_airfoil_geometry move_and_smooth

generate_airfoil_geometry: generate_airfoil_geometry.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -lopenblas -o $@

generate_airfoil_geometry.o: generate_airfoil_geometry.cpp airfoil_bezier.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

move_and_smooth: move_and_smooth.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -lAFEPack -lopenblas -o $@

move_and_smooth.o: move_and_smooth.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f generate_airfoil_geometry generate_airfoil_geometry.o
	rm -f move_and_smooth move_and_smooth.o
