set (TARGETNAME nstream_mexfiles)
set (INCLUDE_DIRS ${nstream_SOURCE_DIR}/src -I${nstream_SOURCE_DIR}/misc/preproc -I${nstream_SOURCE_DIR}/lib/libnstreamclient -I${nstream_SOURCE_DIR}/lib/libutil)
set (DEPENDENCIES ${nstream_SOURCE_DIR}/lib/libnstreamclient/nstreamclient.c ${nstream_SOURCE_DIR}/lib/libnstreamclient/nstreamclient.h ${nstream_SOURCE_DIR}/lib/libutil/exception.cc)
set (LINK_TARGETS ${nstream_BINARY_DIR}/lib/libnstreamclient/libnstreamclient.a ${nstream_BINARY_DIR}/lib/libutil/libutil.a)
set (MEXFILES 
    nstream_getdio
    nstream_getdio_ts
    nstream_getjoy_comedi
    nstream_getmuint_nspike_16
    nstream_getmu_nspike_16
    nstream_getraw_comedi_all
    nstream_getraw_comedi
    nstream_getraw_comedi_decimate
    nstream_getrawint_comedi
    nstream_getrawint_comedi_decimate
    nstream_getrawint_nspike_16
    nstream_getrawint_nspike
    nstream_getrawint_nspike_decimate
    nstream_getraw_nspike_16
    nstream_getraw_nspike_all
    nstream_getraw_nspike_decimate
    nstream_gettime_comedi
    nstream_gettime_nspike
    nstream_exit
    nstream_set_dac_channel
    nstream_set_dac_gain
    nstream_set_filters
    nstream_set_verbose
    nstream_set_channel
    nstream_add_auxdsp
    nstream_add_masterdsp
    nstream_start_record
    nstream_stop_record
    nstream_start_acquire
    nstream_enable_dataglove
    nstream_get_dataglove
    nstream_get_dataglove_ts
    nstream_getmu_nspike_4
    nstream_getmuint_nspike_4
    nstream_getraw_nspike_4
    nstream_getrawint_nspike_4
)

# setup file extension for architecture
#
IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Regular x86
    SET(MEXEXT mexglx)
ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # AMD64:
    SET(MEXEXT mexa64)
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)


# generate a build step for each mex file
#
foreach (MEXFILE ${MEXFILES})
    # add build command
    ADD_CUSTOM_COMMAND(
        OUTPUT    ${CMAKE_CURRENT_BINARY_DIR}/${MEXFILE}.${MEXEXT}
        DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}.c ${DEPENDENCIES}
        COMMAND   mex -cxx
        ARGS      ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}.c ${LINK_TARGETS} -I${INCLUDE_DIRS} -outdir ${CMAKE_CURRENT_BINARY_DIR} 
        COMMENT   ${MEXFILE} 
      )
    # add command to build documentation .m file
    ADD_CUSTOM_COMMAND(
        OUTPUT    ${CMAKE_CURRENT_BINARY_DIR}/${MEXFILE}.m
        DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}.c
        COMMAND   grep ^% ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}.c >${CMAKE_CURRENT_BINARY_DIR}/${MEXFILE}.m
        COMMENT   ${MEXFILE}.m
      )

      set (MEX_BINARIES ${MEX_BINARIES}
            ${CMAKE_CURRENT_BINARY_DIR}/${MEXFILE}.${MEXEXT}
      )

      set (MEX_DOCS ${MEX_DOCS}
            ${CMAKE_CURRENT_BINARY_DIR}/${MEXFILE}.m
      )

endforeach (MEXFILE)

# add build targets for all the built mexfiles 
#
add_custom_target( ${TARGETNAME} ALL
    DEPENDS ${MEX_BINARIES} ${MEX_DOCS}
)

install(FILES ${MEX_BINARIES} ${MEX_DOCS} DESTINATION mex)
