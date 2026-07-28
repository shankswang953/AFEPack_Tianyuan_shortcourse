include $(dir $(lastword $(MAKEFILE_LIST)))../course_config.mk

CXXFLAGS := -O2 -g -std=c++20
LDLIBS := -lAFEPack -lopenblas

.PHONY: all clean

all: main

main: two_heater_goal_adaptive.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

two_heater_goal_adaptive.o: two_heater_goal_adaptive.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f main two_heater_goal_adaptive.o
