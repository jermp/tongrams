cmake_minimum_required(VERSION 2.8)
project(TONGRAMS)

SET(CMAKE_INCLUDE_CURRENT_DIR ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

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

if(TONGRAMS_USE_SANITIZER)
  if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")
endif()

find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_subdirectory(emphf)

add_executable(build_mph_lm build_mph_lm.cpp)
target_link_libraries(build_mph_lm ${Boost_LIBRARIES})

add_executable(build_trie_lm build_trie_lm.cpp)
target_link_libraries(build_trie_lm ${Boost_LIBRARIES})

add_executable(lookup_perf_test test/lookup_perf_test.cpp)
target_link_libraries(lookup_perf_test ${Boost_LIBRARIES})

add_executable(print_stats utils/print_stats.cpp)
target_link_libraries(print_stats ${Boost_LIBRARIES})

add_executable(sort_grams sorters/sort_grams.cpp)
target_link_libraries(sort_grams ${Boost_LIBRARIES})

add_executable(sort_arpa sorters/sort_arpa.cpp)
target_link_libraries(sort_arpa ${Boost_LIBRARIES})

add_executable(check_count_model test/check_count_model.cpp)
target_link_libraries(check_count_model ${Boost_LIBRARIES})

add_executable(score score.cpp)
target_link_libraries(score ${Boost_LIBRARIES})

add_executable(hash_compact_vector_test test/hash_compact_vector_test.cpp)
target_link_libraries(hash_compact_vector_test ${Boost_LIBRARIES})

add_executable(compact_vector_test test/compact_vector_test.cpp)
target_link_libraries(compact_vector_test ${Boost_LIBRARIES})

add_executable(fast_ef_sequence_test test/fast_ef_sequence_test.cpp)
target_link_libraries(fast_ef_sequence_test ${Boost_LIBRARIES})

add_executable(uniform_pef_sequence_test test/uniform_pef_sequence_test.cpp)
target_link_libraries(uniform_pef_sequence_test ${Boost_LIBRARIES})
