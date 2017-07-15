
#TARGET_SOURCES(
#    hipgisaxs
#    PRIVATE
SET(ffcpu_SOURCES
	${CMAKE_CURRENT_LIST_DIR}/ff.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_ana.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_num.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_box.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_cube.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_cylinder.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ff_ana_prism.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ff_ana_prism3x.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ff_ana_prism6.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_pyramid.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_sphere.cpp
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_cone.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ff_ana_sawtooth.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cpu/ff_num_cpu.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cpu/ff_tri_cpu.cpp
)

ADD_LIBRARY(ffcpu ${ffcpu_SOURCES})
