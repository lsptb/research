add_subdirectory(ext)
#add_subdirectory(tests)

########### next target ###############

SET(nspike_LIB_SRCS
   netkeepalive.cc
   nspike.cc
   dsp.cc
   masterdsp.cc
   udpsocket.cc
   auxdsp.cc
   channel.cc
   daqengine.cc
   ringbuffer.cc
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/libcfu-0.03/include ${nstream_SOURCE_DIR}/lib/libutil)

add_library(nspike STATIC ${nspike_LIB_SRCS})

target_link_libraries(nspike rt m cfu util)

set_target_properties(nspike PROPERTIES VERSION 1.0.0 SOVERSION 1)
#install(TARGETS nspike DESTINATION lib/libnspike)


########### install files ###############




#original Makefile.am contents follow:

### Process this file with automake to produce Makefile.in
#
#lib_LTLIBRARIES = libnspike.la
#libnspike_la_SOURCES = nspike.c dsp.c masterdsp.c udpsocket.c auxdsp.c channel.c daqengine.c exception.c ringbuffer.c netkeepalive.c
#libnspike_la_LDFLAGS = -pthread -lrt -lm -static
#libnspike_la_LIBADD = ../ext/libcfu-0.03/src/libcfu.a 
#
#libnspike_la_CFLAGS = -I../ext/libcfu-0.03/include -O0
