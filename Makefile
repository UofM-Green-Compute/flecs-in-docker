#!/bin/bash

# Base Working Directory
BWD := $(shell pwd)

BWDMOUNT := -v $(BWD):$(BWD):ro
BUILDMOUNT := -v $(BWD)/build:$(BWD)/build:delegated
BINMOUNT := -v $(BWD)/bin:$(BWD)/bin:delegated

MOUNTS := $(BWDMOUNT) $(BUILDMOUNT) $(BINMOUNT)

all:
	@echo "make docker - build docker container"
	@echo "make prepare - create build location"
	@echo "make dockerbash - run bash inside the container"
	@echo "make dockerbuild - build the code inside the container"
	@echo "make clean - wipe the build"
	@echo
	@echo "NB: final artefacts live in 'bin'"

docker:
	docker build -t buildenv -f Dockerfile .

prepare:
	mkdir -p $(BWD)/build

clean:
	rm -rf $(BWD)/build
	rm -rf $(BWD)/bin/*

dockerbash: prepare
	docker run -it --rm -u1000:1000 -e BWD=$(BWD) \
	           $(BWDMOUNT) \
	           $(BUILDMOUNT) \
	           $(BINMOUNT) \
	           buildenv \
	           /bin/bash


dockerbuild: prepare
	docker run -it --rm -u1000:1000 -e BWD=$(BWD) \
	           $(MOUNTS) \
	           buildenv \
	           make -f $(BWD)/src/Makefile

devloop:
	make clean
	make prepare
	make dockerbuild
	make BWD=$(BWD) -f $(BWD)/src/Makefile run
