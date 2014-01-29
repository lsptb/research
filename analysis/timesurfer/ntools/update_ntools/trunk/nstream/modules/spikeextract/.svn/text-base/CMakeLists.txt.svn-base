
add_subdirectory(mex)
########### next target ###############

SET(module_spikeextract_SRCS
   module_spikeextract.c
)

add_executable(module_spikeextract ${module_spikeextract_SRCS})

target_link_libraries(module_spikeextract)

install(TARGETS module_spikeextract DESTINATION bin)


########### install files ###############




#original Makefile.am contents follow:

### Process this file with automake to produce Makefile.in
#
#bin_PROGRAMS = module_spikeextract
#module_spikeextract_SOURCES = module_spikeextract.c $(top_builddir)/nstream/nstream.h
#module_spikeextract_CFLAGS  = $(MATLAB_FLAGS)
#module_spikeextract_LDFLAGS = $(WRAPS) $(MATLAB_LINK)
#module_spikeextract_LDADD   = $(MATLAB_LIB)
#
