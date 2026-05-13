find_package(GTest REQUIRED)

set(CMAKE_TEST_OUTPUT_DIRECTORY Tests_SaiModel)

# add unit tests
add_executable(tests
    Tests_SaiModel.cpp
)

target_link_libraries(tests
    PRIVATE ${SAI-MODEL_LIBRARIES}
    PRIVATE ${SAI-URDF_LIBRARIES}
	PRIVATE GTest::GTest
	PRIVATE GTest::Main
)

set_target_properties(tests PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
)

# copy the tests urdfs
add_custom_command(
  TARGET tests
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  "${CMAKE_CURRENT_SOURCE_DIR}/urdf"
  "${CMAKE_CURRENT_BINARY_DIR}/urdf"
)

add_custom_command(
  TARGET tests
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  "${PROJECT_SOURCE_DIR}/urdf_models/HRP4C"
  "${CMAKE_CURRENT_BINARY_DIR}/HRP4C"
)

# Enable testing with CTest
enable_testing()

# Add a test that runs the test executable
add_test(NAME Tests_SaiModel
    COMMAND tests
)
