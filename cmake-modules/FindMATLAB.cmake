# - Locate MATLAB executables
#
# Variables set by this module:
#
#  MATLAB_FOUND
#  MATLAB_EXECUTABLE
#  MATLAB_MEX_EXECUTABLE
#  MATLAB_MEXEXT

include(FindPackageHandleStandardArgs)

find_program(MATLAB_EXECUTABLE matlab DOC "MATLAB program")
mark_as_advanced(MATLAB_EXECUTABLE)

list(FIND MATLAB_FIND_COMPONENTS mex _MEX_REQUESTED)
if(_MEX_REQUESTED GREATER -1)
  if(MATLAB_EXECUTABLE)
    get_filename_component(_MATLAB_REALPATH ${MATLAB_EXECUTABLE} REALPATH)
    get_filename_component(_MATLAB_BINDIR ${_MATLAB_REALPATH} PATH)

    find_program(MATLAB_MEX_EXECUTABLE mex NAMES mex.bat
                 HINTS ${_MATLAB_BINDIR}
                 NO_DEFAULT_PATH
                 DOC "MEX compiler")
    mark_as_advanced(MATLAB_MEX_EXECUTABLE)

    find_program(_MEXEXT_COMMAND mexext NAMES mexext.bat
                 HINTS ${_MATLAB_BINDIR})
    if(_MEXEXT_COMMAND)
      execute_process(OUTPUT_VARIABLE MATLAB_MEXEXT
                      COMMAND ${_MEXEXT_COMMAND}
                      OUTPUT_STRIP_TRAILING_WHITESPACE)
      set(MATLAB_MEXEXT ${MATLAB_MEXEXT} CACHE STRING
          "MEX file extension" FORCE)
      mark_as_advanced(MATLAB_MEXEXT)
    endif()
    unset(_MEXEXT_COMMAND CACHE)
  endif()

  find_package_handle_standard_args(MATLAB DEFAULT_MSG
                                    MATLAB_EXECUTABLE
                                    MATLAB_MEX_EXECUTABLE
                                    MATLAB_MEXEXT)
else()
  find_package_handle_standard_args(MATLAB DEFAULT_MSG MATLAB_EXECUTABLE)
endif()

function(_mex_compile SOURCE BUILDDIR FLAGS OBJ_VAR)
  if(NOT IS_ABSOLUTE ${SOURCE})
    set(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE})
  endif()
  file(RELATIVE_PATH SOURCE_REL ${CMAKE_CURRENT_SOURCE_DIR} ${SOURCE})
  get_filename_component(OUTDIR ${SOURCE_REL} PATH)
  string(REPLACE ".." "__" OUTDIR "${OUTDIR}")
  if(OUTDIR)
    set(OUTDIR ${BUILDDIR}/${OUTDIR})
  else()
    set(OUTDIR ${BUILDDIR})
  endif()
  file(MAKE_DIRECTORY ${OUTDIR})
  get_filename_component(BASENAME ${SOURCE} NAME_WE)
  if(CMAKE_HOST_WIN32)
    set(OBJ ${OUTDIR}/${BASENAME}.obj)
  else()
    set(OBJ ${OUTDIR}/${BASENAME}.o)
  endif()
  file(RELATIVE_PATH OBJ_REL ${CMAKE_BINARY_DIR} ${OBJ})
  set(${OBJ_VAR} ${OBJ} PARENT_SCOPE)

  
  add_custom_command(OUTPUT ${OBJ}
    COMMAND ${MATLAB_MEX_EXECUTABLE} ${SOURCE} -largeArrayDims -c
      -outdir ${OUTDIR} ${FLAGS}
    DEPENDS ${SOURCE} IMPLICIT_DEPENDS ${SOURCE} VERBATIM
    COMMENT "Building MEX object ${OBJ_REL}")
endfunction()

function(add_mex MEXFILE)
  if(NOT MATLAB_MEX_EXECUTABLE)
    message(SEND_ERROR "mex command not found")
  endif()

  set(SOURCES)
  set(FLAGS)
  set(active_list SOURCES)
  foreach(arg ${ARGN})
    if(arg STREQUAL "FLAGS")
      set(active_list FLAGS)
    else()
      list(APPEND ${active_list} "${arg}")
    endif()
  endforeach()

  set(SKIP_OPTIMIZED NO)
  set(SKIP_DEBUG YES)
  if(CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo|MinSizeRel")
    set(FLAGS -O ${FLAGS})
  endif()
  if(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo|Debug")
    set(FLAGS -g ${FLAGS})
    set(SKIP_OPTIMIZED YES)
    set(SKIP_DEBUG NO)
  endif()

  set(BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${MEXFILE}.dir)
  set(OBJECTS)
  set(SKIP NO)
  foreach(arg ${SOURCES})
    if(NOT SKIP)
      if(arg STREQUAL "optimized")
        set(SKIP ${SKIP_OPTIMIZED})
      elseif(arg STREQUAL "debug")
        set(SKIP ${SKIP_DEBUG})
      elseif(arg MATCHES "\\.c$|\\.cpp$|\\.cc$|\\.cxx$|\\.C$")
        _mex_compile(${arg} ${BUILDDIR} "${FLAGS}" OBJ)
        list(APPEND OBJECTS ${OBJ})
      else()
        list(APPEND OBJECTS ${arg})
      endif()
    else()
      set(SKIP NO)
    endif()
  endforeach()

  set(MEXFILE_EXT ${MEXFILE}.${MATLAB_MEXEXT})
  if(NOT CMAKE_HOST_WIN32)
    set(CXX_FLAG "-cxx")
    set(OBJECTS_ARG ${OBJECTS})
  else()
    set(RSPFILE ${BUILDDIR}/objects-arg)
    set(action WRITE)
    foreach(obj ${OBJECTS})
      file(${action} ${RSPFILE} "${obj}\n")
      set(action APPEND)
    endforeach()
    set(OBJECTS_ARG @${RSPFILE})
  endif()
  if(SKIP_DEBUG AND CMAKE_STRIP)
    if(CMAKE_HOST_APPLE)
      set(STRIP_COMMAND COMMAND ${CMAKE_STRIP} -x ${MEXFILE_EXT})
    else()
      set(STRIP_COMMAND COMMAND ${CMAKE_STRIP} --strip-unneeded ${MEXFILE_EXT})
    endif()
  endif()
  add_custom_command(OUTPUT ${MEXFILE_EXT}
    COMMAND ${MATLAB_MEX_EXECUTABLE} -largeArrayDims ${CXX_FLAG}
      -output ${MEXFILE} -outdir ${CMAKE_CURRENT_BINARY_DIR}
      ${OBJECTS_ARG} ${FLAGS}
    ${STRIP_COMMAND}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${OBJECTS} VERBATIM
    COMMENT "Linking MEX library ${MEXFILE_EXT}")
  add_custom_target(${MEXFILE} ALL DEPENDS ${MEXFILE_EXT})
endfunction()
