
SET(sfcpu_SOURCES
	${CMAKE_CURRENT_LIST_DIR}/liquid.cpp
	${CMAKE_CURRENT_LIST_DIR}/sf.cpp 
)

ADD_LIBRARY(sfcpu ${sfcpu_SOURCES})
