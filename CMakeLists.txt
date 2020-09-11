cmake_minimum_required(VERSION 3.15)
project(projectAgNew C)

set(CMAKE_C_STANDARD 90)

add_executable(projectAgNew main.c utils.c utils.h spmat.c spmat.h divider.c divider.h)