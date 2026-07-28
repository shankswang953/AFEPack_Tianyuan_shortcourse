include $(dir $(lastword $(MAKEFILE_LIST)))../course_config.mk

CXXFLAGS := -O2 -g -std=c++20
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
