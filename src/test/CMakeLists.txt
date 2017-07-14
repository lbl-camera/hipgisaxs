
TARGET_SOURCES( 
    hipgisaxs
    PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/qgrid_test_create.cpp
	${CMAKE_CURRENT_LIST_DIR}/test_conv.cpp
	${CMAKE_CURRENT_LIST_DIR}/test_ff.cpp
	${CMAKE_CURRENT_LIST_DIR}/test_image.cpp
	${CMAKE_CURRENT_LIST_DIR}/test_read.cpp
)
