# Shared compile/link configuration for the AFEPack teaching examples.
#
# Direct `make` calls and all run.sh launchers use these values.  A local
# course_config.local file is optional and is shared with course_config.sh.

AFEPACK_EXAMPLES_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
-include $(AFEPACK_EXAMPLES_ROOT)/course_config.local

CXX ?= c++

# AFEPack: finite-element assembly, solves, refinement, and mesh operations.
AFEPACK_PREFIX ?= $(HOME)
OPENBLAS_PREFIX ?= /opt/local
BOOST_INCLUDE ?= /opt/local/libexec/boost/1.81/include

AFEPACK_CPPFLAGS := -I$(AFEPACK_PREFIX)/include \
                    -I$(OPENBLAS_PREFIX)/include/openblas \
                    -I$(BOOST_INCLUDE)
AFEPACK_LDFLAGS := -L$(AFEPACK_PREFIX)/lib \
                   -L$(OPENBLAS_PREFIX)/lib \
                   -Wl,-rpath,$(OPENBLAS_PREFIX)/lib

CPPFLAGS += $(AFEPACK_CPPFLAGS)
LDFLAGS += $(AFEPACK_LDFLAGS)
