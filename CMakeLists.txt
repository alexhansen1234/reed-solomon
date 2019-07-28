project(reed_solomon LANGUAGES CXX)
cmake_minimum_required(VERSION 3.12)

if( NOT ${CMAKE_BINARY_DIR} STREQUAL ${PROJECT_SOURCE_DIR}/build )
  message(FATAL_ERROR "build in the build directory")
endif()

set(CMAKE_CXX_COMPILER /usr/local/bin/gcc-9)

add_executable(${PROJECT_NAME} src/reed_solomon.cpp)

target_include_directories(${PROJECT_NAME}
  PUBLIC
  include
)

target_link_libraries(${PROJECT_NAME}
  PUBLIC
  stdc++
)

target_compile_options(${PROJECT_NAME}
  PUBLIC
  -std=c++2a
  -Wall
  -O3
)

file(GLOB TESTS tests/*)
foreach(INFILE ${TESTS})
  file(RELATIVE_PATH INFILE ${PROJECT_SOURCE_DIR} ${INFILE})
  configure_file(${INFILE} ${INFILE} COPYONLY)
endforeach()
