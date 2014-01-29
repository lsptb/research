
########### next target ###############

SET(util_LIB_SRCS
   exception.cc
)

#include_directories(${nstream_SOURCE_DIR}/src)

add_library(util STATIC ${util_LIB_SRCS})

#target_link_libraries(nspike rt m cfu)

set_target_properties(util PROPERTIES VERSION 1.0.0 SOVERSION 1 COMPILE_FLAGS -fPIC)
#install(TARGETS nstreamclient DESTINATION lib/libnstreamclient)


########### install files ###############

