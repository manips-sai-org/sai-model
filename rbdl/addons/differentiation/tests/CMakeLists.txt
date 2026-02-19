PROJECT (RBDL_ADDON_DIFFERENTIATION_TESTS)

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# Needed for UnitTest++
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../CMake )

# Look for unittest++
FIND_PACKAGE (UnitTest++ REQUIRED)
INCLUDE_DIRECTORIES (${UNITTEST++_INCLUDE_DIR})

SET ( TESTS_SRCS
    runtests.cc
		SpatialAlgebraTest.cc
		rbdl_mathutilsADvsFDTests.cc
		rbdl_utilsADvsFDTests.cc
		DynamicsADvsFDTests.cc
		KinematicsADvsFDTests.cc
    ConstraintsADvsFDTests.cc
    ModelCheckADvsFD.cc
)

INCLUDE_DIRECTORIES ( ../src/ )

SET_TARGET_PROPERTIES ( ${PROJECT_EXECUTABLES} PROPERTIES
  LINKER_LANGUAGE CXX
)

ADD_EXECUTABLE ( rbdl_ad_tests ${TESTS_SRCS} )
# target_compile_features(rbdl_ad_tests PRIVATE cxx_range_for)
SET_TARGET_PROPERTIES ( rbdl_ad_tests PROPERTIES
    LINKER_LANGUAGE CXX
    OUTPUT_NAME runtests
    )

SET (RBDL_LIBRARY rbdl)
IF (RBDL_BUILD_STATIC)
    SET (RBDL_LIBRARY rbdl-static)
ENDIF (RBDL_BUILD_STATIC)

SET (RBDL_LIBRARY_AD rbdl_ad)
IF (RBDL_BUILD_STATIC)
    SET (RBDL_LIBRARY_AD rbdl_ad-static)
ENDIF (RBDL_BUILD_STATIC)

TARGET_LINK_LIBRARIES ( rbdl_ad_tests
        ${UNITTEST++_LIBRARY}
        ${RBDL_LIBRARY}
        ${RBDL_LIBRARY_AD}
    )

IF (RUN_AUTOMATIC_TESTS)
ADD_CUSTOM_COMMAND (TARGET runtests
    POST_BUILD
    COMMAND ./runtests
    COMMENT "Running automated tests..."
    )
ENDIF (RUN_AUTOMATIC_TESTS)
