# executable will have the same name as its directory
get_filename_component(EXEC_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# write the resulting executable in the _current_ directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# add executable from source files in current directory assuming
# to add more source files, manually register them using `add_executable`
get_filename_component(MAIN ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_executable(${EXEC_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/${MAIN}.cpp" glucoseGrid.cpp glucoseTransport.cpp)
find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(${EXEC_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})

# link the executable to `hemocell` and `HDF5` dependencies
target_link_libraries(${EXEC_NAME} ${PROJECT_NAME})
target_link_libraries(${EXEC_NAME} ${HDF5_C_HL_LIBRARIES} ${HDF5_LIBRARIES})

# Standalone numerical tests: no MPI initialization or GoogleTest download.
add_executable(vessel_glucose_test EXCLUDE_FROM_ALL tests/test_glucose_grid.cpp glucoseGrid.cpp)
target_include_directories(vessel_glucose_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(vessel_glucose_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
if(BUILD_TESTING)
    add_test(NAME vessel_glucose_grid COMMAND vessel_glucose_test)
endif()
