# Building the Codebase

This documentation is based on notes targetted at MPhys students, and
assuming very little (if any) understanding of git, Linux, or development
tools.  If that description matches you and you're trying this and it
doesn't work for you, please file an issue with a question.


## Overview

The build is containerised using docker. Practically speaking this means two
things:

* The build is reproducible - no "this works for me on this operating system
  but not you on yours". It builds/breaks for everyone or no-one.

* This means that in order to build the codebase you essentially need 2
  things installed: `make` and `docker`

You need to be set up on github as well - to be able to check out and check
in.

The easiest way of getting these installed on the major platforms is:

* Linux - `apt install build-essential` followed by installing docker from
  get.docker.com

* Mac - install the xcode terminal tools; install docker desktop (but don't
  use docker desktop)

* Windows - use WSL 2 (Windows Subsystem for Linux); then follow Linux
  instructions.


## Clone the repository

I'm assuming you haven't got any development directories or similar set up
yet. If you have, just skip to clone the repository.

Change to your home directory:

    cd ~

Create a directory for development:

    mkdir Development

Move into it:

    cd Development

Clone the repository locally:

    git clone git@github.com:UofM-Green-Compute/flecs-in-docker.git  # Check out the source

Change into that directory

    cd flecs-in-docker

And have a look around:

    ls
    find . |grep -v git 
    tree                 # won't work on mac os out the box

If that all looks good, you're good to go!


## Building / Etc

Once the dependencies are installed, the build process should all "just
work".  Continuing on from being in the above directory, you should now be
able to:

Build the docker container:

    make docker


Build the code inside the docker container:

    make dockerbuild


Run the code inside the docker container:

    make run

Clean up your dev environment:

    make clean

Automatically do "make clean; make devbuild; make run":

    make devloop

Day to day you probably want to edit, make chanes, build and run, so this:

    make dockerbuild
    make run

... is probably what you want.

Doing this build results in an application called `bin/ecs_application` . If
you are using linux or WSL under windows, you can run this using either

    make run

Or this:

    ./bin/ecs_application

Under Mac, you must use `make run` (at least for now).





