# Simple, reproducible Makefile for C++20/23 + Flecs (single-file C lib)
# Works inside Ubuntu 24.04 LTS container with build-essential installed.

APP_BINARY := ecs_application

# Discover base working dir (repo root) from this Makefile’s location
ifndef BWD
	BWD := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/..)
endif

SRC := $(BWD)/src
INC := $(BWD)/include
OBJ := $(BWD)/bin
RUNDIR := $(BWD)/outputs

# --- toolchain & flags -------------------------------------------------------
CXX      ?= g++
CC       ?= gcc
CPPFLAGS  = -I$(INC) -MMD -MP
CXXFLAGS  = -std=gnu++23 -O2 -g -Wall -Wextra -Wpedantic
# CFLAGS    = -std=c11     -O2 -g -Wall -Wextra -Wpedantic -D_POSIX_C_SOURCE=199309L
CFLAGS    = -std=gnu99     -O2 -g -Wall -Wextra -Wpedantic
LDFLAGS   =
LDLIBS    = -pthread

# --- sources & objects -------------------------------------------------------
CXX_SOURCES := $(wildcard $(SRC)/*.cpp)
C_SOURCES   := $(SRC)/flecs.c
CXX_OBJECTS := $(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(CXX_SOURCES))
C_OBJECTS   := $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(C_SOURCES))
OBJECTS     := $(C_OBJECTS) $(CXX_OBJECTS)
DEPS        := $(OBJECTS:.o=.d)

app := $(OBJ)/$(APP_BINARY)

# --- rules -------------------------------------------------------------------
.PHONY: all clean run dirs
all: dirs $(app)

dirs:
	@mkdir -p $(OBJ)

$(app): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

# C++ source
$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# C source (flecs)
$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) -f $(OBJECTS) $(DEPS) $(app)

run: all
	cd $(RUNDIR) ; $(app)

-include $(DEPS)

