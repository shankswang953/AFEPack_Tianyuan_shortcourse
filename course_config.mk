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

# On macOS, compile for the architecture provided by libAFEPack.  This matters
# on Apple Silicon when an existing AFEPack installation was built under
# Rosetta.  A machine-specific COURSE_TARGET_ARCH can still override the
# automatic choice in course_config.local.
ifeq ($(shell uname -s),Darwin)
  AFEPACK_DYLIB := $(firstword \
    $(wildcard $(AFEPACK_PREFIX)/lib/libAFEPack.dylib) \
    $(wildcard $(AFEPACK_PREFIX)/lib/libAFEPack.0.dylib))
  ifeq ($(origin COURSE_TARGET_ARCH),undefined)
    ifneq ($(strip $(AFEPACK_DYLIB)),)
      AFEPACK_DYLIB_INFO := $(shell file "$(AFEPACK_DYLIB)" 2>/dev/null)
      ifeq ($(findstring arm64,$(AFEPACK_DYLIB_INFO)),)
        ifneq ($(findstring x86_64,$(AFEPACK_DYLIB_INFO)),)
          COURSE_TARGET_ARCH := x86_64
        endif
      endif
    endif
  endif
  ifneq ($(strip $(COURSE_TARGET_ARCH)),)
    COURSE_ARCH_FLAGS := -arch $(COURSE_TARGET_ARCH)
  endif
endif

AFEPACK_CPPFLAGS := -I$(AFEPACK_PREFIX)/include \
                    -I$(OPENBLAS_PREFIX)/include/openblas \
                    -I$(BOOST_INCLUDE)
AFEPACK_LDFLAGS := -L$(AFEPACK_PREFIX)/lib \
                   -L$(OPENBLAS_PREFIX)/lib \
                   -Wl,-rpath,$(OPENBLAS_PREFIX)/lib

CPPFLAGS += $(AFEPACK_CPPFLAGS) $(COURSE_ARCH_FLAGS)
LDFLAGS += $(AFEPACK_LDFLAGS) $(COURSE_ARCH_FLAGS)
