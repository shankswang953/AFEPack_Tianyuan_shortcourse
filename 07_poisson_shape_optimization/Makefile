include $(dir $(lastword $(MAKEFILE_LIST)))../course_config.mk

CXXFLAGS := -O2 -g -std=c++20 -Wall -Wextra
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
