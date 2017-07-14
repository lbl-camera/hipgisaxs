
TARGET_SOURCES(
    hipgisaxs
    PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/compute_params.cpp
	${CMAKE_CURRENT_LIST_DIR}/inst_detector.cpp
	${CMAKE_CURRENT_LIST_DIR}/inst_scattering.cpp
	${CMAKE_CURRENT_LIST_DIR}/layer.cpp
	${CMAKE_CURRENT_LIST_DIR}/multilayer.cpp
	${CMAKE_CURRENT_LIST_DIR}/qgrid.cpp
	${CMAKE_CURRENT_LIST_DIR}/shape.cpp
	${CMAKE_CURRENT_LIST_DIR}/structure.cpp
	${CMAKE_CURRENT_LIST_DIR}/unitcell.cpp
)
