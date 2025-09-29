#pragma once

#include <coroutine>
#include <stdexcept>


template<typename T>
class simple_generator {
public:
    struct promise_type;  
    using handle_type = std::coroutine_handle<promise_type>;

private:
    handle_type mCoro;

public:
    explicit simple_generator(handle_type h) : mCoro(h) {}

    simple_generator(simple_generator &&other_sg) noexcept : mCoro(other_sg.mCoro) {
        other_sg.mCoro = nullptr;
    }
    simple_generator &operator=(simple_generator &&other) noexcept {
        if (this != other) {
            mCoro = other.mCoro;
            other.mCoro = nullptr;
            return *this;
        }
    }
    simple_generator(const simple_generator &) = delete;
    simple_generator &operator=(const simple_generator &) = delete;
    ~simple_generator() { if (mCoro) { mCoro.destroy(); } }

    // Implementation of the external API called by the user to actually use the generator
    void start() {  try_next(); }
    bool running() { 
        return not mCoro.done();
    }
    void try_next() { 
        mCoro.resume();
        if (mCoro.promise().m_latest_exception) {
            std::rethrow_exception(mCoro.promise().m_latest_exception);
        }
    }
    T take() { return std::move(mCoro.promise().m_current_value); }

    // Implementation of the internal API called when co_yield/etc are triggered inside the coroutine
    class promise_type {
        T m_current_value;
        std::exception_ptr m_latest_exception;
        friend simple_generator;

    public:
        auto get_return_object()      { return simple_generator{ handle_type::from_promise(*this) } ; }
        auto yield_value(T some_value) {
            m_current_value = some_value;   // Capture the yielded value
            return std::suspend_always{};
        }
        auto unhandled_exception() { m_latest_exception = std::current_exception(); }
        auto initial_suspend()        { return std::suspend_always{}; }
        auto final_suspend() noexcept { return std::suspend_always{}; }
         auto return_void()           { return std::suspend_never{}; }
    };

private:
    // Implementation of the iterator protocol
    class iterator {
        simple_generator<T> &owner;
        bool done;
        void iter_next() {
            owner.try_next();
            done = not owner.running();
        }
    public:
        bool operator != (const iterator &r) const { return done != r.done; }
        auto operator * () const { return owner.take(); }
        iterator &operator ++ () { iter_next(); return *this; }
        iterator(simple_generator<T> &o, bool d) : owner(o), done(d) { if ( not done ) iter_next(); }

    };

public:
    // Public access to the internal iterator protocol

    iterator begin() { return iterator{*this, false}; }
    iterator end() { return iterator{*this, true}; }

};
