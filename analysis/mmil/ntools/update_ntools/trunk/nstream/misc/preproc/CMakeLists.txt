
SET(nspike_preproc_SRCS
   nspike_preproc.c
)

SET(comedi_preproc_SRCS
   comedi_preproc.c
)

SET(nspike_preproc_128_SRCS
   nspike_preproc_128.c
)

SET(nspike_preproc_16_SRCS
   nspike_preproc_16.c
)

SET(nspike_preproc_4_SRCS
   nspike_preproc_4.c
)

add_executable(nspike_preproc ${nspike_preproc_SRCS})
add_executable(comedi_preproc ${comedi_preproc_SRCS})
add_executable(nspike_preproc_128 ${nspike_preproc_128_SRCS})
add_executable(nspike_preproc_16 ${nspike_preproc_16_SRCS})
add_executable(nspike_preproc_4 ${nspike_preproc_4_SRCS})

#target_link_libraries(nstream nspike cfu comedi)

install(TARGETS nspike_preproc
                comedi_preproc
                nspike_preproc_128
                nspike_preproc_16
                nspike_preproc_4
        DESTINATION bin)


