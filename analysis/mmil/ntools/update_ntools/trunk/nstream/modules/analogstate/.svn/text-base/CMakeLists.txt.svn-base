add_subdirectory(mex)


########### next target ###############

SET(module_analogstate_SRCS
   module_analogstate.c
)

include_directories(${NSTREAM_SOURCE_DIR}/nstream)

add_executable(module_analogstate ${module_analogstate_SRCS})

target_link_libraries(module_analogstate)

install(TARGETS module_analogstate DESTINATION bin)


########### install files ###############




#original Makefile.am contents follow:

### Process this file with automake to produce Makefile.in
#
#bin_PROGRAMS = module_analogstate
#module_analogstate_SOURCES = module_analogstate.c $(top_builddir)/nstream/nstream.h
#
#SUBDIRS = mex
