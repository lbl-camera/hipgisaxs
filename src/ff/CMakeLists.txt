
TARGET_SOURCES(
    hipgisaxs
    PRIVATE
	src/ff/ff.cpp
	src/ff/ff_ana.cpp
	src/ff/ff_num.cpp
	src/ff/ff_ana_box.cpp
	src/ff/ff_ana_cube.cpp
	src/ff/ff_ana_cylinder.cpp
    src/ff/ff_ana_prism.cpp
    src/ff/ff_ana_prism3x.cpp
    src/ff/ff_ana_prism6.cpp
	src/ff/ff_ana_pyramid.cpp
	src/ff/ff_ana_sphere.cpp
	src/ff/ff_ana_cone.cpp
    src/ff/ff_ana_sawtooth.cpp
    src/ff/cpu/ff_num_cpu.cpp
    src/ff/cpu/ff_tri_cpu.cpp
)

IF (CUDA_FOUND)
    ADD_SUBDIRECTORY(gpu)
    ADD_DEFINITIONS(-DFF_ANA_GPU -DFF_NUM_GPU)
ENDIF (CUDA_FOUND)
