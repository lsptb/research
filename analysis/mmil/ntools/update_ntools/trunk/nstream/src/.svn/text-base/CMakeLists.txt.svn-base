include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib/libnspike ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libnspike/ext/libcfu-0.03/include ${nstream_SOURCE_DIR}/lib/libutil ${nstream_SOURCE_DIR}/lib/libnstreamclient)

add_subdirectory(mex)

########### next target ###############

SET(nstream_SRCS
   nstream.cc
   nstream.h
   dataglovedaq.cc
   dataglovedaq.h
)

add_executable(nstream ${nstream_SRCS})

target_link_libraries(nstream nspike cfu comedi fglove)

install(TARGETS nstream DESTINATION bin)


########### install files ###############


#original Makefile.am contents follow:

### Process this file with automake to produce Makefile.in
#
#bin_PROGRAMS = nstream 
#
#nstream_SOURCES = nstream.c nstream.h $(top_builddir)/src/libnspike.h
#nstream_LDADD = $(top_builddir)/lib/libnspike/libnspike.la $(top_builddir)/lib/libnspike/ext/libcfu-0.03/src/libcfu.a 
#nstream_LDFLAGS = -pthread -lcomedi
#nstream_CFLAGS = -I../lib/libnspike/ext/libcfu-0.03/include
#
