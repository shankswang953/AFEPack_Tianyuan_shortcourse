include $(dir $(lastword $(MAKEFILE_LIST)))../course_config.mk

CXXFLAGS := -O2 -g -std=c++20
LDLIBS := -lAFEPack -lopenblas

.PHONY: all clean

all: main

main: boundary_flux_dual.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

boundary_flux_dual.o: boundary_flux_dual.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f main boundary_flux_dual.o
