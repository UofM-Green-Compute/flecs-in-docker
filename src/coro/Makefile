APP_BINARY := printer

# Note BWD is injected in via docker/host environment - build inside docker due to g++ versions needed
# If unset, being run outside the container - which may not work
# due to dev environment

# nb BWD == Base Working Directory

ifndef BWD
	# This gets the path for the makefile, gets the parent directories, and canonicalises it
	BWD := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/.. )
endif

# Loosely based on:
# https://stackoverflow.com/questions/25966411/autodependencies-for-make-generate-the-d-files-but-are-not-reading-them

BIN := $(BWD)/bin
SRC := $(BWD)/src
INC := $(BWD)/include
OBJ := $(BWD)/build

app     := $(BIN)/$(APP_BINARY)
sources := $(wildcard $(SRC)/*.cpp)
objects := $(subst $(SRC),$(OBJ),$(sources:.cpp=.o))
deps    := $(objects:.o=.d)

CXX      := g++
CPPFLAGS := -I $(INC) -MMD -MP
CXXFLAGS := -std=c++20  -fcoroutines
LDFLAGS  := -L /usr/lib/x86_64-linux-gnu -static
#LDLIBS   := -lcurl
LDLIBS   := 


$(app) : $(objects)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean: ; $(RM) -f $(objects) $(deps) $(app)

-include $(deps)

run:
	$(app)
