cmake_minimum_required( VERSION 2.8.8 )
project(data_loader CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#
set(SPADES_MAIN_SRC_DIR ${CMAKE_SOURCE_DIR}/../../spades/src/)
set(SPADES_MAIN_INCLUDE_DIR ${SPADES_MAIN_SRC_DIR}/include)
set(SPADES_MAIN_DEBRUIJN_DIR ${SPADES_MAIN_SRC_DIR}/debruijn)
set(SPADES_MAIN_ONLINE_VIS_DIR ${SPADES_MAIN_SRC_DIR}/online_vis)
#
set(SPADES_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../spades/build/release)
set(SPADES_TOOLS_BINARY_DIR ${SPADES_BINARY_DIR}/bin)
set(SPADES_BIN ${CMAKE_SOURCE_DIR}/../../spades/bin)
set(SPADES_BUILT_INCLUDE_DIR ${SPADES_BINARY_DIR}/include)
#
include_directories(${SPADES_BUILT_INCLUDE_DIR})
#
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${SPADES_MAIN_DEBRUIJN_DIR})
include_directories(${SPADES_MAIN_ONLINE_VIS_DIR})

include_directories(${SPADES_BIN})
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
#
configure_file("${SPADES_MAIN_INCLUDE_DIR}/config.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/config.hpp")
#
##include_directories(${EXT_DIR}/include/ConsensusCore)
##target_link_libraries(spades)
#
add_executable(data_loader data_loader.cpp )

#target_link_libraries(data_loader spades)
#add_library(spades STATIC IMPORTED)
#SET_PROPERTY(TARGET spades PROPERTY IMPORTED_LOCATION ${SPADES_BIN}/spades)
find_package (SPAdes REQUIRED)
if (SPAdes_FOUND)
  include_directories(${SPAdes_INCLUDE_DIRS})
  target_link_libraries (data_loader ${SPAdes_LIBRARIES})
endif (SPAdes_FOUND)
