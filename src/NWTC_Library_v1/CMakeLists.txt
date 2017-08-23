
set(NWTCLIBS_SOURCES NWTC_Library.f90 NWTC_IO.f90 NWTC_Num.f90 NWTC_Aero.f90
  ModMesh.f90)

if (APPLE OR UNIX)
  if (FCNAME MATCHES "gfortran.*")
    list (APPEND NWTCLIBS_SOURCES SysGnuLinux.f90)
  elseif (FCNAME MATCHES "ifort.*")
    list (APPEND NWTCLIBS_SOURCES SysIFL.f90)
  endif (FCNAME MATCHES "gfortran.*")
elseif (WIN32)
  if (FCNAME MATCHES "gfortran.*")
    list (APPEND NWTCLIBS_SOURCES SysGnuWin.f90)
  elseif (FCNAME MATCHES "ifort.*")
    list (APPEND NWTCLIBS_SOURCES SysIVF.f90)
  endif (FCNAME MATCHES "gfortran.*")
endif (APPLE OR UNIX)

if (DOUBLE_PRECISION)
  list (APPEND NWTCLIBS_SOURCES DoubPrec.f90)
else (DOUBLE_PRECISION)
  list (APPEND NWTCLIBS_SOURCES SingPrec.f90)
endif (DOUBLE_PRECISION)

add_library(nwtclibs  ${NWTCLIBS_SOURCES})
target_link_libraries(nwtclibs ${LAPACK_LIBRARIES} ${CMAKE_DL_LIBS})

install(TARGETS nwtclibs
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
