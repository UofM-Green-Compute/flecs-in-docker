
#include <iostream>
#include <coro.hpp>
#include <futureprint.hpp>

// Couple of simple generators

simple_generator<int> myrange(int max) {
    print("simple_generator<int>\n");
    for(int i=0; i< max; i++) {
        print("simple_generator<for> {} \n", i);
        co_yield i;
    }
}

simple_generator<int> fibs(int max) {
    print("fib generator\n");
    int a {1};
    int b {1};
    for(int i=0; i< max; i++) {
        int n;
         print("fib generator (for) : {}  {} \n", i, max);
        co_yield a;
        n = a + b;
        a = b;
        b = n;
    }
}
