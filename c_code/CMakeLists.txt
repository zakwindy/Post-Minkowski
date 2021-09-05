cmake_minimum_required(VERSION 3.0)

# set the project name
project(BHB VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)

add_executable(binary
	src/main.cpp
	include/boost
	)

add_executable(make_structs
	src2/main2.cpp
	)
