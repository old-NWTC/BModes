
cmake_minimum_required(VERSION 2.8.12)
project(BMODES Fortran)

# Set default installation directory
get_filename_component(default_install
  "${CMAKE_SOURCE_DIR}/../install" ABSOLUTE)
set(CMAKE_INSTALL_PREFIX ${default_install} CACHE string
  "Choose installation directory; default bmodes/install" FORCE)
option(BUILD_SHARED_LIBS "Enabled building shared libraries" off)

# Set default build type to be RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the build type: Debug Release MinSizeRel" FORCE)
endif (NOT CMAKE_BUILD_TYPE)


# Set the fortran modules to be stored in a common directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/ftnmods )
# Add that to the include path
include_directories( ${CMAKE_Fortran_MODULE_DIRECTORY} )

get_filename_component(FCNAME ${CMAKE_Fortran_COMPILER} NAME)

# By default we set the preprocessor directive
if (FCNAME MATCHES "gfortran.*")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
elseif(FCNAME MATCHES "ifort.*")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
endif (FCNAME MATCHES "gfortran.*")

# Set the position independent code directive on architectures that support it
if (APPLE OR UNIX)
   if (FCNAME MATCHES "gfortran.*")	# bjj: this probably works on intel fortran for mac/linux, but I have not tested
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
   endif (FCNAME MATCHES "gfortran.*")
endif (APPLE OR UNIX)

# Force free-form compilation
if (FCNAME MATCHES "gfortran.*")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif (FCNAME MATCHES "gfortran.*")


# use static libraries
if (WIN32)
  if (FCNAME MATCHES "ifort.*")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /libs:static")
  endif (FCNAME MATCHES "ifort.*")
   
endif(WIN32)
   


add_subdirectory(NWTC_Library_v1)
add_subdirectory(bmodes)

include(CPack)

