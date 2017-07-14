# Locate yaml-cpp
#
# This module defines
# yaml-cpp_FOUND, if false, do not try to link to yaml-cpp
# yaml-cpp_LIBS, where to find yaml-cpp
# yaml-cpp_INCLUDE_DIR, where to find yaml.h
#
# By default, the dynamic libraries of yaml-cpp will be found. To find the static ones instead,
# you must set the yaml-cpp_STATIC_LIBRARY variable to TRUE before calling find_package(YamlCpp ...).
#
# If yaml-cpp is not installed in a standard path, you can use the yaml-cpp_DIR CMake variable
# to tell CMake where yaml-cpp is.

include (CheckIncludeFileCXX)
include (CheckCXXSourceRuns)

# find the yaml-cpp include directory
find_path(yaml-cpp_INCLUDE_DIR yaml-cpp/yaml.h
    PATH_SUFFIXES include
    PATHS
    ~/Library/Frameworks/yaml-cpp/include/
    /Library/Frameworks/yaml-cpp/include/
    /usr/local/include/yaml-cpp/
    /usr/local/include/
    /usr/include/yaml-cpp/
    /usr/include/
    /sw/yaml-cpp/ # Fink
    /opt/local/yaml-cpp/ # DarwinPorts
    /opt/csw/yaml-cpp/ # Blastwave
    /opt/yaml-cpp/
    ${yaml-cpp_DIR}/include/)

set(CMAKE_REQUIRED_INCLUDES ${yaml-cpp_INCLUDE_DIR})
set(CMAKE_REQUIRED_QUIET True)

# first look for outdated yaml-cpp0.3 include files
unset(yaml-cpp_FOUND_03 CACHE)
check_include_file_cxx("yaml-cpp/aliasmanager.h" yaml-cpp_FOUND_03)
if(${yaml-cpp_FOUND_03})
    message(WARNING "Found include file for libyaml-cpp0.3. Most probably this precludes libyaml-cpp0.5 from being properly installed")
endif()

# now look for needed yaml-cpp0.5 include files
unset(yaml-cpp_FOUND_05 CACHE)
check_include_file_cxx("yaml-cpp/node/detail/iterator_fwd.h" yaml-cpp_FOUND_05)
if(${yaml-cpp_FOUND_05})
else()
    message(FATAL_ERROR "Include file for libyaml-cpp0.5 not found")
endif()

# attempt to find static library first if this is set
if(yaml-cpp_STATIC_LIBRARY)
    set(yaml-cpp_STATIC libyaml-cpp.a)
endif()

# find the yaml-cpp library
find_library(yaml-cpp_LIBRARY
    NAMES ${yaml-cpp_STATIC} yaml-cpp
    PATH_SUFFIXES lib64 lib
    PATHS ~/Library/Frameworks
    /Library/Frameworks
    /usr/local
    /usr
    /sw
    /opt/local
    /opt/csw
    /opt
    ${yaml-cpp_DIR}/lib)

# try to compile, link, and run a test program
unset(yaml-cpp_RUNS CACHE)
set(CMAKE_REQUIRED_LIBRARIES yaml-cpp)
check_cxx_source_runs("#include \"yaml-cpp/yaml.h\"\n#include <assert.h>\nint main() {\n    YAML::Node node = YAML::Load(\"[1, 2, 3]\");\n    assert(node.IsSequence());\n}" yaml-cpp_RUNS)
if(${yaml-cpp_RUNS})
else()
    message(FATAL_ERROR "Test of libyaml-cpp0.5 failed")
endif()

# handle the QUIETLY and REQUIRED arguments and set yaml-cpp_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(yaml-cpp DEFAULT_MSG yaml-cpp_INCLUDE_DIR yaml-cpp_LIBRARY)
mark_as_advanced(yaml-cpp_INCLUDE_DIR yaml-cpp_LIBRARY)
set(yaml-cpp_LIBRARYS ${yaml-cpp_LIBRARY})
set(yaml-cpp_INCLUDE_DIRS ${yaml-cpp_INCLUDE_DIR})
