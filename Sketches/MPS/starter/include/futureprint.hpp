#pragma once

#include <format>
#include <string>


std::string fmt(const std::string_view str_fmt, auto&&... args) {
    return std::vformat(str_fmt, std::make_format_args(args...));
}

// Simple implementation of a C++23 style print command. Not quite the same,
// but close enough to be useful

// Next line requires g++ 13.x or higher
constexpr void print(const std::string_view str_fmt, auto&&... args) {
    fputs(std::vformat(str_fmt, std::make_format_args(args...)).c_str(),
        stdout);
}
