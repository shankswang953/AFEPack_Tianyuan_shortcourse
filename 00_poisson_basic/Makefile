include $(dir $(lastword $(MAKEFILE_LIST)))../course_config.mk

CXXFLAGS := -O2 -g -std=c++20
LDLIBS := -lAFEPack -lopenblas

.PHONY: all clean

all: main

main: poisson_equation.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

poisson_equation.o: poisson_equation.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f main poisson_equation.o
