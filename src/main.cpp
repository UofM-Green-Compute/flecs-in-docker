#include <coro.hpp>
#include <generator.hpp>
#include <futureprint.hpp>

int main(int, char**) {

    auto ranger = myrange(10);

    for(ranger.start(); ranger.running(); ranger.try_next() ) {
        int i = ranger.take();
        print("Explicit looping returned: {}\n", i);
    }

    for(auto i : myrange(5) ) {
        print("For value in myrange: {}\n", i);
    }

    for(auto i : fibs(7) ) {
        print("FIB : {}\n", i);
    }

    return 0;
}
