set (TARGETNAME misc_mexfiles)
set (INCLUDE_DIRS ${nstream_SOURCE_DIR}/src)
set (DEPENDENCIES ${nstream_SOURCE_DIR}/src/nstream.h)
set (MEXFILES 
    openserial
    closeserial
    readserial
    writeserial
    setparallel
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
        COMMAND   mex 
        ARGS      ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}.c -I${INCLUDE_DIRS} -outdir ${CMAKE_CURRENT_BINARY_DIR} 
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

