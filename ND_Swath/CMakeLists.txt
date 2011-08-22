cmake_minimum_required(VERSION 2.4)

SET(APP_NAME ItkPath-ND_Swath)

PROJECT ( ${APP_NAME} )

# Look for ITK
INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (USE_ITK_FILE)
  INCLUDE(${USE_ITK_FILE})
ENDIF(USE_ITK_FILE)

LINK_LIBRARIES(
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  ITKMetaIO
)

SET(TEST_SRC_DIR ./)

SET(TEST_SRCS ${TEST_SRC_DIR}itkSwathChainCodePathFilterTest.cxx)

ADD_EXECUTABLE( ${APP_NAME}  main.cxx ${TEST_SRCS} )
TARGET_LINK_LIBRARIES( ${APP_NAME} ${LINK_LIBRARIES} )

#SET(TEST_NAME itkChainCodePath2DTest)
#ADD_EXECUTABLE( ${TEST_NAME}  ${TEST_NAME}.cxx)
#TARGET_LINK_LIBRARIES( ${TEST_NAME} ${LINK_LIBRARIES} )

IF(UNIX)
  # -ggdb will CRASH gcc w/ an internal compiler error
  ADD_DEFINITIONS(-Wall -g)
ENDIF(UNIX)
