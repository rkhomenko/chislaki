cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(chislaki VERSION 3.5.0 LANGUAGES CXX)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
   "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -isystem")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DDEBUG")
else()
    # Realization of MSVC compiler support
    message(FATAL_ERROR "You are seriously, man? Use GNU or Clang compiler or go away!")
endif()

set(HEADER_SUFFIX "hpp")
set(SOURCE_SUFFIX "cpp")
set(INCLUDE_DIR "include")
set(SOURCE_DIR "src")
set(TEST_DIR "tests")
set(t)

file(GLOB_RECURSE INCLUDES "${INCLUDE_DIR}/*.${HEADER_SUFFIX}")

include_directories(${INCLUDE_DIR})

# Tests
if(MAKE_TESTS)
    file(GLOB_RECURSE TEST_FILES "${TEST_DIR}/*.${SOURCE_SUFFIX}")
    foreach(TEST_FILE ${TEST_FILES})
        get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME)
        string(REGEX REPLACE ".${SOURCE_SUFFIX}" "" TARGET_NAME ${TEST_FILE_NAME})

        add_executable(${TARGET_NAME} ${TEST_FILE} ${INCLUDES})
        set_property(TARGET ${TARGET_NAME} PROPERTY CXX_STANDARD 17)
    endforeach()
endif()
