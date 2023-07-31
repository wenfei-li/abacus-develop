list(APPEND objects
    LCAO_dftu_new.cpp
    LCAO_dftu_new_fdelta.cpp
    LCAO_dftu_new_gedm.cpp
    LCAO_dftu_new_io.cpp
    LCAO_dftu_new_mpi.cpp
    LCAO_dftu_new_pdm.cpp
    LCAO_dftu_new_psialpha.cpp
    LCAO_dftu_new_vdelta.cpp
    LCAO_dftu_new_interface.cpp
)

add_library(
    dftu_new
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(dftu_new)
endif()

#if(BUILD_TESTING)
#    add_subdirectory(test)
#endif()

