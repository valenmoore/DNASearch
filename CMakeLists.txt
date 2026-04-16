cmake_minimum_required(VERSION 4.1)
project(DNASearch)

set(CMAKE_CXX_STANDARD 20)

find_package(cpr CONFIG REQUIRED)

add_executable(DNASearch main.cpp
        FastaParser.cpp
        FastaParser.h
        Search.cpp
        Search.h)

target_link_libraries(DNASearch PRIVATE cpr::cpr)