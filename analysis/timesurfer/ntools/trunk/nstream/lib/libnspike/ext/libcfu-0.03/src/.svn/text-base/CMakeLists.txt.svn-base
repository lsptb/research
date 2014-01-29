include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include )


########### next target ###############

SET(libcfu_SRCS
cfu.c
cfuconf.c
cfuhash.c
cfulist.c
cfuopt.c
cfustring.c
cfuthread_queue.c
cfutime.c
)
add_library(cfu ${libcfu_SRCS})
SET_TARGET_PROPERTIES( cfu PROPERTIES COMPILE_FLAGS -fPIC)

target_link_libraries(cfu)

#install(TARGETS cfu DESTINATION bin)


########### install files ###############

#install(FILES  $(top_srcdir)/include/cfu.h $(top_srcdir)/include/cfuhash.h      $(top_srcdir)/include/cfutime.h $(top_srcdir)/include/cfustring.h      $(top_srcdir)/include/cfulist.h $(top_srcdir)/include/cfuconf.h      $(top_srcdir)/include/cfuthread_queue.h $(top_srcdir)/include/cfutypes.h      $(top_srcdir)/include/cfuopt.h DESTINATION /include)



#original Makefile.am contents follow:

#lib_LIBRARIES = libcfu.a
#libcfu_a_SOURCES = cfuhash.c cfutime.c cfustring.c cfulist.c cfuconf.c cfuthread_queue.c cfu.c cfuopt.c
#bin_PROGRAMS = libcfu-config
#
#CFU_CONFIG_DEFINES = -DCFU_LIBDIR="\"@libdir@\"" -DCFU_LIBS="\"@LIBS@\"" -DCFU_INCLUDEDIR="\"@includedir@\"" -DCFU_PREFIX="\"@prefix@\""  -DCFU_EXEC_PREFIX="\"@exec_prefix@\""
#
#if DEBUG
#AM_CFLAGS = -Wall -W -DCFU_DEBUG=1 $(CFU_CONFIG_DEFINES)
#else
#AM_CFLAGS = -Wall -W -fPIC $(CFU_CONFIG_DEFINES)
#endif
#
#include_HEADERS = $(top_srcdir)/include/cfu.h $(top_srcdir)/include/cfuhash.h \
#    $(top_srcdir)/include/cfutime.h $(top_srcdir)/include/cfustring.h \
#    $(top_srcdir)/include/cfulist.h $(top_srcdir)/include/cfuconf.h \
#    $(top_srcdir)/include/cfuthread_queue.h $(top_srcdir)/include/cfutypes.h \
#    $(top_srcdir)/include/cfuopt.h
#
## libcfu_a_LDFLAGS = -lpthread
#INCLUDES = -I$(top_srcdir)/include
#
#default: libcfu.a libcfu-config
#
#libcfu-config.o: libcfu-config.c
#
#libcfu-config: libcfu.a libcfu-config.o
#	$(CC) $(CFLAGS) -o libcfu-config libcfu-config.o libcfu.a -pthread
#
