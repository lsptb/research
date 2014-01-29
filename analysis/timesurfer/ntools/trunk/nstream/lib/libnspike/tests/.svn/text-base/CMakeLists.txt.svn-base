
########### next target ###############
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../ext/libcfu-0.03/include ${CMAKE_CURRENT_SOURCE_DIR}/../ ${nstream_SOURCE_DIR}/lib/libutil)

SET(check_libnspike_SRCS
   check_libnspike.c
   check_nspike.c
   check_masterdsp.c
   check_auxdsp.c
   check_channel.c
   check_daqengine.c
   check_ringbuffer.c
   fakehw.c
)

add_executable(check_libnspike ${check_libnspike_SRCS})

target_link_libraries(check_libnspike nspike cfu check)


########### install files ###############




#original Makefile.am contents follow:

### Process this file with automake to produce Makefile.in
#
#TESTS = check_libnspike
#check_PROGRAMS = check_libnspike
#check_libnspike_SOURCES = check_libnspike.c check_nspike.c check_masterdsp.c check_auxdsp.c check_channel.c check_daqengine.c check_ringbuffer.c fakehw.c $(top_builddir)/src/libnspike.h
#check_libnspike_CFLAGS = @CHECK_CFLAGS@ -I../ext/libcfu-0.03/include
#check_libnspike_LDADD = @CHECK_LIBS@ $(top_builddir)/src/libnspike.la ../ext/libcfu-0.03/src/libcfu.a 
#check_libnspike_LDFLAGS = -pthread
