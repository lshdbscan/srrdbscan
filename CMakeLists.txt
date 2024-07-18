cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(LSHDBSCAN)

if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-std=c++17 -O3  -mavx2")
#set(CMAKE_CXX_FLAGS "-std=c++17 -O3 -ffast-math -mavx2")
message(STATUS "**************Build type: ${CMAKE_BUILD_TYPE}************")

if(NOT DEFINED HIGHFIVE_USE_BOOST)
  set(HIGHFIVE_USE_BOOST Off)
endif()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

add_subdirectory(third_party/HighFive)


include_directories(src/include)
file (GLOB LIB_SOURCES src/*.cc)
file (GLOB LIB_HEADERS src/include/*.h)
add_library(myLib ${LIB_SOURCES} ${LIB_HEADERS})
target_link_libraries(myLib HighFive)

add_executable(SRR_LSHDBSCAN_exec apps/main_SRR_DBSCAN.cc)
#target_link_libraries(SRR_LSHDBSCAN_exec myLib gsl gslcblas gomp tbb pthread)
target_link_libraries(SRR_LSHDBSCAN_exec myLib gomp tbb pthread)


# for testing that HighFive works
#add_executable(writeTest apps/fileWritingTest.cc)
#target_link_libraries(writeTest HighFive)


