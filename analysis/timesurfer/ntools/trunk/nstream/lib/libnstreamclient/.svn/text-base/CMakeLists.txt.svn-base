
########### next target ###############

SET(nstreamclient_LIB_SRCS
   nstreamclient.c
)

include_directories(${nstream_SOURCE_DIR}/src ${nstream_SOURCE_DIR}/lib/libutil)

add_library(nstreamclient STATIC ${nstreamclient_LIB_SRCS})

#target_link_libraries(nspike rt m cfu)
set_target_properties(nstreamclient PROPERTIES VERSION 1.0.0 SOVERSION 1 COMPILE_FLAGS -fPIC)
#install(TARGETS nstreamclient DESTINATION lib/libnstreamclient)


########### install files ###############

