include $(dir $(lastword $(MAKEFILE_LIST)))../course_config.mk

CXXFLAGS := -O2 -g -std=c++20
LDLIBS := -lAFEPack -lopenblas

.PHONY: all clean

all: tianyuan_mesh_merge mesh_merge_demo poisson_solver

tianyuan_mesh_merge: tianyuan_mesh_merge.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

tianyuan_mesh_merge.o: tianyuan_mesh_merge.cpp local_mesh_merge.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

mesh_merge_demo: mesh_merge_demo.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

mesh_merge_demo.o: mesh_merge_demo.cpp local_mesh_merge.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

poisson_solver: poisson_equation.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

poisson_equation.o: poisson_equation.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f tianyuan_mesh_merge tianyuan_mesh_merge.o \
	      mesh_merge_demo mesh_merge_demo.o \
	      poisson_solver poisson_equation.o
