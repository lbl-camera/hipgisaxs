/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num.hpp
 *  Created: Nov 05, 2011
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _NUMERIC_FF_HPP_
#define _NUMERIC_FF_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <cstring>

#ifdef USE_MPI
#include <mpi.h>
#include <woo/comm/multi_node_comm.hpp>
#endif

#include <common/parameters.hpp>
#include <common/enums.hpp>
#ifdef USE_HDF5
#include <file/hdf5shape_reader.h>
#endif
#include <model/qgrid.hpp>
#include <utils/utilities.hpp>
#ifdef FF_NUM_GPU
  #include <ff/gpu/ff_num_gpu.cuh>  // for gpu version
#elif defined USE_MIC
  #include <ff/mic/ff_num_mic.hpp>  // for mic version
#else
  #include <ff/cpu/ff_num_cpu.hpp>  // for cpu version
#endif // FF_NUM_GPU

#include <numerics/matrix.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>

namespace hig {

  /**
   * Some declarations.
   */
  #ifdef FF_NUM_GPU
    void write_slice_to_file(cucomplex_t *ff, int nqx, int nqy, int nqz,
                  char* filename, int axis, int slice);
  #endif
  
  
  /**
   * Class for computing Form Factor across multiple nodes using MPI.
   */
  class NumericFormFactor {
        // TODO: make the cpu and cpu version classes children of this class ...
    public:

      #ifdef FF_NUM_GPU    // use GPUs for numerical
        #ifdef FF_NUM_GPU_FUSED
          NumericFormFactor(int block_cuda_y, int block_cuda_z):
                  block_cuda_t_(0),
                  block_cuda_y_(block_cuda_y), block_cuda_z_(block_cuda_z),
                  gff_(block_cuda_y, block_cuda_z) { }
        #endif
        #ifdef KERNEL2
          NumericFormFactor(int block_cuda_t, int block_cuda_y, int block_cuda_z):
                  block_cuda_t_(block_cuda_t), block_cuda_y_(block_cuda_y),
                  block_cuda_z_(block_cuda_z),
                  gff_(block_cuda_t, block_cuda_y, block_cuda_z) { }
        #else
          NumericFormFactor(int block_cuda): block_cuda_(block_cuda), gff_(block_cuda) { }
        #endif // KERNEL2
      #elif defined USE_MIC  // use MICs for numerical
        NumericFormFactor(): mff_() { }
      #else          // use CPUs for numerical
        NumericFormFactor(): cff_() { }
      #endif  // FF_NUM_GPU

      ~NumericFormFactor() { }

      bool init(RotMatrix_t &, std::vector<complex_t>&);
      void clear() { }        // TODO ...

      bool compute(const char* filename, std::vector<complex_t>& ff,
              RotMatrix_t &
              #ifdef USE_MPI
                , woo::MultiNode&, std::string
              #endif
              );
  
      bool compute2(const char* filename, std::vector<complex_t>& ff,
              RotMatrix_t &
              #ifdef USE_MPI
                , woo::MultiNode&, std::string
              #endif
              );
    private:

      // TODO: make these for gpu only ...
      unsigned int block_cuda_;
      /* cuda block size for main kernel */
      unsigned int block_cuda_t_;
      unsigned int block_cuda_y_;
      unsigned int block_cuda_z_;

      #ifdef FF_NUM_GPU
        NumericFormFactorG gff_;    // for computation on the GPU
      #elif defined USE_MIC
        NumericFormFactorM mff_;    // for computation on the MIC and/or intel CPU
      #else
        NumericFormFactorC cff_;    // for computation only on CPU
      #endif

      unsigned int nqx_;
      unsigned int nqy_;
      unsigned int nqz_;

      RotMatrix_t rot_;
  
      ShapeFileType get_shapes_file_format(const char*);
      unsigned int read_shapes_file_dat(const char* filename, real_vec_t& shape_def);
      unsigned int read_shapes_file(const char* filename,
//                    #ifndef __SSE3__
                      real_vec_t& shape_def
//                    #else
//                      #ifdef USE_GPU
//                        real_vec_t &shape_def
//                      #else
//                        real_t* &shape_def
//                      #endif
//                    #endif
                    );
      void find_axes_orientation(std::vector<real_t> &shape_def, std::vector<short int> &axes);
      bool construct_ff(int p_nqx, int p_nqy, int p_nqz,
                int nqx, int nqy, int nqz,
                int p_y, int p_z,
                #ifdef FF_NUM_GPU
                  cucomplex_t* p_ff,
                #else
                  complex_t* p_ff,
                #endif
                complex_vec_t& ff,
                #ifdef USE_MPI
                  woo::MultiNode&, std::string,
                #endif
                real_t&, real_t&);

  }; // class NumericFormFactor
} // namespace hig

#endif // _NUMERIC_FF_HPP_
