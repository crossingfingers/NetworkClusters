cmake_minimum_required(VERSION 3.15)
project(cProject C)

set(CMAKE_C_STANDARD 90)

add_executable(cProject main.c spmat.c divider.c stack.c stack.h)