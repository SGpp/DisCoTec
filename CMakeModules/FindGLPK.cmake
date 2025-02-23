find_path(
  GLPK_INCLUDE_DIR
  glpk.h
)

if( GLPK_INCLUDE_DIR )
  find_library(
    GLPK_LIBRARY
    NAMES glpk
  )
  if( GLPK_LIBRARY )
    set(GLPK_LIBRARY_DIR "")
    get_filename_component(GLPK_LIBRARY_DIRS ${GLPK_LIBRARY} PATH)
    set(GLPK_FOUND ON)
    set(GLPK_INCLUDE_DIRS ${GLPK_INCLUDE_DIR})
    set(GLPK_LIBRARIES ${GLPK_LIBRARY})
  endif(GLPK_LIBRARY)
else(GLPK_INCLUDE_DIR)
  message(STATUS "FindGLPK: Could not find glpk.h")
endif(GLPK_INCLUDE_DIR)

if(GLPK_FOUND)
  if(NOT GLPK_FIND_QUIETLY)
    message(STATUS "FindGLPK: Found both glpk.h and libglpk.a")
  endif(NOT GLPK_FIND_QUIETLY)
else(GLPK_FOUND)
  if(GLPK_FIND_REQUIRED)
    message(FATAL_ERROR "FindGLPK: Could not find glpk.h and/or libglpk.a")
  endif(GLPK_FIND_REQUIRED)
endif(GLPK_FOUND)
