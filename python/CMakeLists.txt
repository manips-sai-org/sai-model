set(PYBIND11_FINDPYTHON ON)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

add_library(sai_model_py SHARED bindings.cpp)
target_compile_features(sai_model_py PRIVATE cxx_std_17)
target_link_libraries(sai_model_py PRIVATE pybind11::headers Python::Python)
target_link_options(sai_model_py PRIVATE "-Wl,-flat_namespace")

target_link_libraries(sai_model_py PRIVATE ${SAI-MODEL_LIBRARIES}
	${SAI-URDF_LIBRARIES})
target_include_directories(sai_model_py PRIVATE
	${SAI-MODEL_INCLUDE_DIRS}
	${SAI-URDF_INCLUDE_DIRS}
	${PARSER_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIR}
	${RBDL_INCLUDE_DIR}
	${ED_RBDL_INCLUDE_DIR}
)

set_target_properties(sai_model_py PROPERTIES
	PREFIX ""
	SUFFIX ".${Python_SOABI}.so"
)
