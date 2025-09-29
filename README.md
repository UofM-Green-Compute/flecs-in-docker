This docker repo was designed initially for use with C++20 coroutines and
various other tests. It's being repurposed as a build environment that
abstracts away the flecs dependencies.



Old notes:

Despite appearances, this directory is not actually about docker. It's about
using docker to test the availability of C++20 formats in g++ in current versions
of ubuntu.
---------

Log:

dockerbuild - updated to noble (24.04) rather than mantic (23.10)

- Reason is that when this sketch was created 23.10 was a future release
  relative to the last LTS (22.04LTS). 22.04's gcc was incomplete with
  regard to C++23 support. I didn't want to update the host from a current
  long term support version (LTS), so built the container instead - which
  can run future versions. However, 23.10 is now less recent than noble, so
  have simply updated to the latest LTS instead.

