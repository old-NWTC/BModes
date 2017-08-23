
add_executable(bmodes Alloc.f90 BldVib.f90 BModes.f90
  eigsolver.f heap_sort.f90 modules.f90 SetProg.f90 Struct.f90 utilss.f90)
target_link_libraries(bmodes nwtclibs)

install(TARGETS bmodes
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib )
