cmake_minimum_required(VERSION 3.0)

# set the project name
project(binary VERSION 1.0)

set(SOURCES main.cpp body.cpp include/body.h)

add_executable(binary SOURCES)
