cmake_minimum_required(VERSION 2.8)
set(TARG_PFX mgt-glm-)

set(COMMON_SRCS delcher.cc fasta.cc gene.cc)
set(ICM_SRCS icm.cc)
set(ICM_PROGS build-icm build-fixed score-fixed simple-score)

foreach(TARG_NAME ${ICM_PROGS})
    
    set(TARG_NAME_PREF ${TARG_PFX}${TARG_NAME})
    add_executable(${TARG_NAME_PREF}  ${TARG_NAME}.cc) 
    target_link_libraries(${TARG_NAME_PREF} GLMicm GLMcommon m)
    list(APPEND ALL_EXECS ${TARG_NAME_PREF})

endforeach(TARG_NAME)


add_library(GLMcommon ${COMMON_SRCS})
add_library(GLMicm ${ICM_SRCS})

INSTALL(TARGETS ${ALL_EXECS} RUNTIME DESTINATION bin)

#SET(CMAKE_BUILD_TYPE distribtion)
#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")

#SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
#SET(CMAKE_C_FLAGS_DEBUG -g "-O0")
