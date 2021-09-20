cmake_minimum_required(VERSION 2.8)
project(TONGRAMS)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  MESSAGE(STATUS "DEBUG defined")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
endif()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

if(TONGRAMS_USE_POPCNT)
  if(UNIX)
    # Use popcount intrinsic. Available on x86-64 since SSE4.2.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTONGRAMS_USE_POPCNT")
  endif()
endif()

if(TONGRAMS_USE_PDEP)
  if(UNIX)
    # Use parallel bit deposit intrinsic.
    # Available on x86-64 since Intel’s Haswell CPUs: BMI2 (Bit Manipulation Instruction Set 2).
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTONGRAMS_USE_PDEP")
  endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")

  if(TONGRAMS_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

find_package(Boost COMPONENTS iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${TONGRAMS_SOURCE_DIR}/include)

add_subdirectory(external/emphf)

file(GLOB SRC_SOURCES src/*.cpp)
foreach(SRC ${SRC_SOURCES})
  get_filename_component (SRC_NAME ${SRC} NAME_WE) # without extension
  add_executable(${SRC_NAME} ${SRC})
  target_link_libraries(${SRC_NAME} ${Boost_LIBRARIES})
endforeach(SRC)

file(GLOB TEST_SOURCES test/test_*.cpp)
foreach(TEST_SRC ${TEST_SOURCES})
  get_filename_component (TEST_SRC_NAME ${TEST_SRC} NAME_WE) # without extension
  add_executable(${TEST_SRC_NAME} ${TEST_SRC})
  target_link_libraries(${TEST_SRC_NAME} ${Boost_LIBRARIES})
endforeach(TEST_SRC)
