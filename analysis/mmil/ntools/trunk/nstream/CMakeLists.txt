#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

project(nstream)

cmake_minimum_required(VERSION 2.6)
set (CMAKE_INSTALL_PREFIX ${nstream_SOURCE_DIR}/pkg)

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(modules)
add_subdirectory(misc)

########### install files ###############




#original Makefile.am contents follow:

### Process this file with automake to produce Makefile.in
#
#SUBDIRS = modules 
