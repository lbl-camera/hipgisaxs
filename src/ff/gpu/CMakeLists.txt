
#TARGET_SOURCES(
#    hipgisaxs
#    PRIVATE
SET(ffgpu_SOURCES
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_box_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_cube_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_cylinder_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_prism6_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_prism_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_pyramid_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_ana_sphere_gpu.cu
	${CMAKE_CURRENT_LIST_DIR}/ff_tri_gpu.cu
)

CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
CUDA_ADD_LIBRARY(ffgpu ${ffgpu_SOURCES})
