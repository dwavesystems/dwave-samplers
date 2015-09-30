find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(SWIG REQUIRED)

include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

set(PYTHON_CXXFLAGS)
set(_PYTHON_SYMBOL_FILE ${CMAKE_BINARY_DIR}/python-symbols)
if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
  set(PYTHON_LIBRARIES)
  file(WRITE ${_PYTHON_SYMBOL_FILE} "_init*\n")
  set(PYTHON_CXXFLAGS -fvisibility=hidden)
  set(PYTHON_LDFLAGS
    -Wl,-undefined,dynamic_lookup
    -Wl,-exported_symbols_list,${_PYTHON_SYMBOL_FILE})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  set(PYTHON_LIBRARIES)
  file(WRITE ${_PYTHON_SYMBOL_FILE}
"{
  global: init*;
  local: *;
};")
  set(PYTHON_CXXFLAGS -fvisibility=hidden)
  set(PYTHON_LDFLAGS -Wl,--version-script=${_PYTHON_SYMBOL_FILE})
elseif(MSVC_VERSION AND NOT MSVC_VERSION LESS 1800)
  set(PYTHON_CXXFLAGS /DHAVE_ROUND)
endif()
string(REPLACE ";" " " PYTHON_CXXFLAGS "${PYTHON_CXXFLAGS}")
string(REPLACE ";" " " PYTHON_LDFLAGS "${PYTHON_LDFLAGS}")

if(WIN32)
  set(PYTHON_EXTENSION_PROPERTIES SUFFIX .pyd)
else()
  set(PYTHON_EXTENSION_PROPERTIES SUFFIX .so) # yes, even on OSX
endif()
