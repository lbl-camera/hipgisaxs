
TARGET_SOURCES(
    hipgisaxs
    PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/edf_reader.cpp
	${CMAKE_CURRENT_LIST_DIR}/objectshape_reader.cpp
	${CMAKE_CURRENT_LIST_DIR}/rawshape_reader.cpp
	${CMAKE_CURRENT_LIST_DIR}/read_oo_input.cpp
)

IF (USE_HDF5)
    TARGET_SOURCES(
        hipgisaxs
        PRIVATE
	    ${CMAKE_CURRENT_LIST_DIR}/hdf5shape_reader.c
    )
ENDIF (USE_HDF5)
