include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include )

########### next target ###############

SET(hash_usage_SRCS
   hash_usage.c
)

add_executable(hash_usage ${hash_usage_SRCS})

target_link_libraries(hash_usage cfu)


########### next target ###############

SET(large_hash_SRCS
   large_hash.c
)

add_executable(large_hash ${large_hash_SRCS})

target_link_libraries(large_hash cfu)


########### next target ###############

SET(strings_SRCS
   strings.c
)

add_executable(strings ${strings_SRCS})

target_link_libraries(strings cfu)


########### next target ###############

SET(conf_example_SRCS
   conf_example.c
)

add_executable(conf_example ${conf_example_SRCS})

target_link_libraries(conf_example cfu)


########### next target ###############

SET(conf_example2_SRCS
   conf_example2.c
)

add_executable(conf_example2 ${conf_example2_SRCS})

target_link_libraries(conf_example2 cfu)


########### next target ###############

SET(opt_example_SRCS
   opt_example.c
)

add_executable(opt_example ${opt_example_SRCS})

target_link_libraries(opt_example cfu)


########### install files ###############




#original Makefile.am contents follow:

#noinst_PROGRAMS = hash_usage large_hash strings conf_example conf_example2 opt_example
#LDADD = -L$(top_srcdir)/src -lcfu
#
#if DEBUG
#AM_CFLAGS = -Wall -W -I$(top_srcdir)/include -DCFU_DEBUG=1
#else
#AM_CFLAGS = -Wall -W -I$(top_srcdir)/include
#endif
#
#AM_LDFLAGS =  -Xlinker -rpath -Xlinker $(top_src_dir)/src
#
