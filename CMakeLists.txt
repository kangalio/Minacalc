cmake_minimum_required(VERSION 3.15)
project(minacalc)

set(CMAKE_CXX_STANDARD 14)

add_executable(minacalc main.cpp minacalc.cpp minacalc.h NoteDataStructures.h smloader.cpp smloader.h)