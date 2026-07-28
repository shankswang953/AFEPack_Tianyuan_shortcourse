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

all: indicator_study local_refine poisson_solver

indicator_study: indicator_study.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

indicator_study.o: indicator_study.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

local_refine: local_refine.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

local_refine.o: local_refine.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

poisson_solver: poisson_equation.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

poisson_equation.o: poisson_equation.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f indicator_study indicator_study.o \
	      local_refine local_refine.o poisson_solver poisson_equation.o
