cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 11) 

# if(NOT CMAKE_BUILD_TYPE)
# set(CMAKE_BUILD_TYPE Release)
# endif()


SET(DEBUG 0)

IF(DEBUG MATCHES 1)
    add_compile_options(-g -ggdb -O0 -std=c++11) #debug
ELSE()
    add_compile_options(-O3 -w -fopenmp -std=c++11) #-ggdb
ENDIF()

# project name
project("HMS_AFC")

# Include directories
find_package(yaml-cpp REQUIRED)  # crea le variabili ${YAML_INCLUDE_DIRS} e ${YAML_CPP_LIBRARIES}
message(STATUS "YAML_INCLUDE_DIRS = " ${YAML_INCLUDE_DIRS})
message(STATUS "YAML_CPP_LIBRARIES = " ${YAML_CPP_LIBRARIES})

# Include directories:: WORHP
FIND_LIBRARY(WORHP_LIBRARY NAMES worhp)
FIND_PATH(WORHP_INCLUDE_DIR NAMES worhp/worhp.h)
IF(NOT WORHP_LIBRARY)
    MESSAGE(FATAL_ERROR "WORHP support was requested, but the library could not be located. The library should be called 'libworhp')")
ENDIF(NOT WORHP_LIBRARY)
IF(NOT WORHP_INCLUDE_DIR)
    MESSAGE(FATAL_ERROR "WORHP support was requested, but worhp/worhp.h not be located.")
ENDIF(NOT WORHP_INCLUDE_DIR)
MESSAGE(STATUS "WORHP library: ${WORHP_LIBRARY}")
MESSAGE(STATUS "WORHP include dir: ${WORHP_INCLUDE_DIR}")

# Paths
set (SRCS_EXE   
    src/main.cpp
    src/ut.cpp
    src/aux_covariance.cpp
    src/C_KeplerArc.cpp                
    src/rendezvousUT_options.cpp
    src/C_prb_RR_HMS.cpp
    src/ode.cpp
    # src/C_prb_rendezvousUT.cpp
    # src/C_prb_RR_MB.cpp
    # src/C_prb_RR_HMS_uDV.cpp
    # src/earthmars.cpp
    # RR/src/prb_RR.cpp
)

set (EXE_INCLUDE_DIR include)
#set(CSPICE_INCLUDE_DIR "/usr/local/include/cspice")  # Set the path to CSPICE include directory
#set(CSPICE_LIBRARY "/usr/local/include/cspice/lib/cspice.a")  # Set the path to the CSPICE library file

# set (EOS_HOME ${PROJECT_SOURCE_DIR}/ext/EOS)
# set (EOS_INCLUDE_DIR ${EOS_HOME}/include)

set(EIGEN_INCLUDE_DIR "/usr/local/include/Eigen")

include_directories("/usr/local/include")

# Include
include_directories(    
    ${EXE_INCUDE_DIR}
    ${EIGEN_INCLUDE_DIR}  
    ${EOS_INCLUDE_DIR}
    ${WORHP_INCLUDE_DIR}
    ${YAML_INCLUDE_DIRS}
    #${CSPICE_INCLUDE_DIR}
)

# add the executable
set (EXE_NAME ${PROJECT_NAME})
add_executable(${EXE_NAME} ${SRCS_EXE})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# link: inserisci le librerie da linkare
set (LIBS  
    #  ${EOS_LIBS}  
    # ${LIBS_SUNDIALS}
    ${YAML_CPP_LIBRARIES}
    ${WORHP_LIBRARY}
    #${CSPICE_LIBRARY}
)

# auto-link
target_link_libraries(${EXE_NAME} ${LIBS})
target_link_libraries(${PROJECT_NAME})





 
 

 

 

 