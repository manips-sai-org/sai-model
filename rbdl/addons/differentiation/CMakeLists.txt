PROJECT (RBDL_ADDON_DIFFERENTIATION)

CMAKE_MINIMUM_REQUIRED (VERSION 3.10)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -Wall -Wno-long-long")

LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake )

# FIND_PACKAGE (UnitTest++ REQUIRED)

INCLUDE_DIRECTORIES (
	${BUILD_DIR}/include/rbdl
	${BUILD_DIR}/include/rbdl
	${BINARY_DIR}/include
#   ${UNITTEST++_INCLUDE_DIR}
)

SET_TARGET_PROPERTIES ( ${PROJECT_EXECUTABLES} PROPERTIES
		LINKER_LANGUAGE CXX
	)

# Perform the proper linking
SET (CMAKE_SKIP_BUILD_RPATH FALSE)
SET (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# define add-on sources
SET ( DIFFERENTIATION_SOURCES
	# ModelAD.cc
	ModelED.cc
	# FdModelEntry.cc
	# ModelEntryFDC.cc
	# ConstraintsED.cc
	# ConstraintsAD.cc
	# ConstraintsFD.cc
	# ConstraintsFDC.cc
	KinematicsED.cc
	# KinematicsAD.cc
	# KinematicsFD.cc
	# KinematicsFDC.cc
	# DynamicsAD.cc
	# DynamicsFD.cc
	# DynamicsFDC.cc
	DynamicsED.cc
	# JointAD.cc
	JointED.cc
	# rbdl_utilsAD.cc
	# rbdl_utilsFD.cc
)

# Iterate over each source file and prepend the current source directory
foreach(src ${RBDL_SOURCES})
    list(APPEND RBDL_SOURCES_WITH_DIR "${BUILD_DIR}/${src}")
endforeach()

IF ( RBDL_BUILD_STATIC )
	ADD_LIBRARY ( rbdl_ad-static STATIC
		${RBDL_SOURCES_WITH_DIR}
		${DIFFERENTIATION_SOURCES}
		)
	# target_compile_features(rbdl_ad-static PRIVATE cxx_range_for)
	SET_TARGET_PROPERTIES ( rbdl_ad-static PROPERTIES PREFIX "lib")
	SET_TARGET_PROPERTIES ( rbdl_ad-static PROPERTIES OUTPUT_NAME "rbdl_ad")

	INSTALL (TARGETS rbdl_ad-static
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	)
ELSE (RBDL_BUILD_STATIC)
	ADD_LIBRARY (
		rbdl_ad SHARED
		${RBDL_SOURCES_WITH_DIR}
		${DIFFERENTIATION_SOURCES}
	)
	# target_compile_features(rbdl_ad PRIVATE cxx_range_for)
	SET_TARGET_PROPERTIES ( rbdl_ad PROPERTIES
		VERSION ${RBDL_VERSION}
		SOVERSION ${RBDL_SO_VERSION}
	)

	INSTALL (TARGETS rbdl_ad
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	)
ENDIF (RBDL_BUILD_STATIC)

INCLUDE_DIRECTORIES (
	./
	./../../src/
	./../../include/
)

# FILE ( GLOB headers
# 	"${CMAKE_CURRENT_SOURCE_DIR}/*.h"
# 	)

SET (headers 
	${CMAKE_CURRENT_SOURCE_DIR}/SpatialAlgebraOperatorsED.h 
	${CMAKE_CURRENT_SOURCE_DIR}/JointED.h 
	${CMAKE_CURRENT_SOURCE_DIR}/KinematicsED.h 
	${CMAKE_CURRENT_SOURCE_DIR}/ModelED.h 
	${CMAKE_CURRENT_SOURCE_DIR}/DynamicsED.h)

INSTALL ( FILES ${headers}
	DESTINATION
	${CMAKE_INSTALL_INCLUDEDIR}/rbdl/addons/differentiation
	)

IF (RBDL_BUILD_TESTS)
	ADD_SUBDIRECTORY ( tests )
ENDIF (RBDL_BUILD_TESTS)
