cmake_minimum_required(VERSION 3.16)
project(hicreate LANGUAGES CXX)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(hicreate
    main.cpp
    config.h
    matrix.cpp
    matrix.h
    reference.cpp
    reference.h
    fragmenter.cpp
    fragmenter.h
    simulator.cpp
    simulator.h
)
target_link_libraries(hicreate PRIVATE Threads::Threads)
