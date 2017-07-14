
TARGET_SOURCES(
    hipgisaxs
    PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/matmul.cpp
	${CMAKE_CURRENT_LIST_DIR}/string_utils.cpp
	${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
)
