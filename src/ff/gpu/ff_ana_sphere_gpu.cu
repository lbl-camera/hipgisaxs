/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_gpu.cu
 *  Created: Oct 16, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <fstream>
#include <complex>
#include <cuComplex.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include <ff/gpu/ff_ana_gpu.cuh>
#include <common/enums.hpp>
#include <common/constants.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>
#include <utils/gpu/cu_utilities.cuh>

namespace hig {

  extern __constant__ real_t transvec_d[3];

  /** Form Factor of Sphere:
   *  R : (real) Radius of the sphere 
   *  q : (complex) q-vector
   *  ff = L * W * H exp(j * qz * H/2) * sinc(qx * L/2) * sinc (qy * W/2) * sinc(qz * H/2)
   */
  __device__  __inline__ cucomplex_t FormFactorSphere(cucomplex_t qx, cucomplex_t qy, cucomplex_t qz, 
          real_t radius){
    cucomplex_t qR    = cuCsqrt(qx * qx + qy * qy + qz * qz) * radius;
    if (cuC_abs(qR) < CUTINY_){
        return make_cuC(REAL_ZERO_, REAL_ZERO_);
    }
    real_t     vol = 4 * PI_ * radius * radius * radius;
    cucomplex_t sincos = cuCsin(qR) - cuCcos(qR) * qR; 
    cucomplex_t expval = cuCexpi(qz * radius);
    cucomplex_t qR3 = qR * qR * qR;
    return (vol * expval * sincos / qR3);
  }

  __device__ __inline__ cucomplex_t CutSphere(cucomplex_t qx, cucomplex_t qy, cucomplex_t qz,
          real_t rz, real_t z){
      cucomplex_t qp = cuCsqrt(qx * qx + qy * qy);
      if (cuC_abs(qp) < CUTINY_){
          return make_cuC(REAL_ZERO_, REAL_ZERO_);
      }
      real_t t1 = 2 * PI_ * rz * rz;
      real_t t2 = rz * qp.x;
      if (t2 < CUTINY_)
          return make_cuC(REAL_ZERO_, REAL_ZERO_);
      cucomplex_t expval = cuCexpi(qz * z);
      cucomplex_t p1 = qp * qz;
      // TODO will produce wrong answers for spheres buried inside layers
      cucomplex_t bess = make_cuC(j1(t2)/t2, REAL_ZERO_);
      return (t1 * bess * expval);
  }


__device__ cucomplex_t integrate(cucomplex_t qx, cucomplex_t qy, cucomplex_t qz,
          real_t radius, real_t height, int N) {
    cucomplex_t out = make_cuC(REAL_ZERO_, REAL_ZERO_);
    cucomplex_t rval = make_cuC(REAL_ZERO_, REAL_ZERO_);
    real_t cen = height - radius;
    real_t dh = height /(real_t) (N-1);
    for (int i = 0; i < N; i++) {
      real_t z = (i*dh) - cen;
      real_t rz = sqrt(radius * radius - z * z);
      if (i == 0 || i == N-1)
        out = out + CutSphere(qx, qy, qz, rz, z);
      else
        out = out + 2.f * CutSphere(qx, qy, qz, rz, z);
    }
    return 0.5f * dh * out;
  }

  __global__ void ff_sphere_kernel (unsigned int nqy, unsigned int nqz,
          real_t * qx, real_t * qy, cucomplex_t * qz, cucomplex_t * ff,
          RotMatrix_t rot,
          int nr, real_t * r, real_t * distr_r, 
          int nh, real_t * h, real_t * distr_h){

    int i_z = blockDim.x * blockIdx.x + threadIdx.x;

    // sub-kernel params
    const int n_quadrature_pts = 512;
    if (i_z < nqz){
      int i_y = i_z % nqy;
      cucomplex_t c_neg_unit = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
      cucomplex_t mqx, mqy, mqz;
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      cucomplex_t temp_ff = make_cuC(REAL_ZERO_, REAL_ZERO_);
      for (int i_r = 0; i_r < nr; i_r++)
        for (int i_h = 0; i_h < nh; i_h++) {
          cucomplex_t t1 = cuCexpi(mqz * (h[i_h] - r[i_r]));
          cucomplex_t fval = integrate(mqx, mqy, mqz, r[i_r], h[i_h], n_quadrature_pts);
          temp_ff = temp_ff + (distr_h[i_h] * distr_r[i_r] * t1 * fval);
        }
      cucomplex_t temp = transvec_d[0] * mqx + transvec_d[1] * mqy + transvec_d[2] * mqz;
      ff[i_z] = temp_ff * cuCexpi(temp);
    }
  }


  __global__ void ff_full_sphere_kernel (unsigned int nqy, unsigned int nqz, 
          real_t * qx, real_t * qy, cucomplex_t * qz, cucomplex_t * ff,
          RotMatrix_t rot,
          int nr, real_t * r, real_t * distr_r){
    int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    if (i_z < nqz){
      int i_y = i_z % nqy;
      cucomplex_t c_neg_unit = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
      cucomplex_t mqx, mqy, mqz;
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      cucomplex_t temp_ff = make_cuC(REAL_ZERO_, REAL_ZERO_);
      for (int i = 0; i < nr; i++) {
        temp_ff = temp_ff + distr_r[i] * FormFactorSphere(mqx, mqy, mqz, r[i]);
      }
      cucomplex_t temp1 = transvec_d[0] * mqx + transvec_d[1] * mqy + transvec_d[2] * mqz;
      ff[i_z] =  temp_ff * cuCexpi(temp1);
    }
  } // ff_full_sphere_kernel()

   
  bool AnalyticFormFactorG::compute_sphere(
                  const std::vector<real_t>& x,
                  const std::vector<real_t>& distr_x,
                  const std::vector<real_t>& h,
                  const std::vector<real_t>& distr_h,
                  const RotMatrix_t & rot, const std::vector<real_t>& transvec,
                  std::vector<complex_t>& ff) {

    unsigned int n_x = x.size(), n_distr_x = distr_x.size();
    unsigned int n_h = h.size(), n_distr_h = distr_h.size();
    const real_t *x_h = x.empty() ? NULL : &*x.begin();
    const real_t *distr_x_h = distr_x.empty() ? NULL : &*distr_x.begin();
    const real_t *h_h = h.empty() ? NULL : &h[0];
    const real_t *distr_h_h = distr_h.empty() ? NULL : &distr_h[0];
    real_t transvec_h[3] = {transvec[0], transvec[1], transvec[2]};

    // construct device buffers
    real_t *x_d, *distr_x_d;
    real_t *h_d, *distr_h_d;

    cudaMalloc((void**) &x_d, n_x * sizeof(real_t));
    cudaMalloc((void**) &distr_x_d, n_distr_x * sizeof(real_t));

    // copy data to device buffers
    cudaMemcpy(x_d, x_h, n_x * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_x_d, distr_x_h, n_distr_x * sizeof(real_t), cudaMemcpyHostToDevice);

    // if the shphere is truncated, copy extra parameters to device memory
    if (n_h > 0) {
      cudaMalloc((void**) &h_d, n_h * sizeof(real_t));
      cudaMalloc((void**) &distr_h_d, n_distr_h * sizeof(real_t));
      cudaMemcpy(h_d, h_h, n_x * sizeof(real_t), cudaMemcpyHostToDevice);
      cudaMemcpy(distr_h_d, distr_h_h, n_distr_h * sizeof(real_t), cudaMemcpyHostToDevice);
    }
    
    cudaMemcpyToSymbol(transvec_d, transvec_h, 3*sizeof(real_t), 0, cudaMemcpyHostToDevice); 

    int num_threads = 256;
    int num_blocks =  nqz_ / num_threads + 1;
    dim3 ff_grid_size(num_blocks, 1, 1);
    dim3 ff_block_size(num_threads, 1, 1);
    //std::cerr << "Q-Grid size = " << nqz_ << std::endl;

    // the kernel
    if (n_h > 0) {
      ff_sphere_kernel <<<num_blocks, num_threads >>> (nqy_, nqz_, 
            qx_, qy_, qz_, ff_, rot,
            n_x, x_d, distr_x_d, n_h, h_d, distr_h_d);
    } else {
      ff_full_sphere_kernel<<<num_blocks, num_threads>>> (nqy_, nqz_, 
            qx_, qy_, qz_, ff_, rot,
            n_x, x_d, distr_x_d);
    }
    
    construct_output_ff(ff);

    /**
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            std::cout << ff[i*100 + j].real() << "    ";
        }
        std::cout << std::endl;
    }
    std::exit(7);
    **/
    if (n_h > 0) {
      cudaFree(distr_h_d);
      cudaFree(h_d);
    }
    cudaFree(distr_x_d);
    cudaFree(x_d);
   
    return true;
  } // AnalyticFormFactorG::compute_sphere()

} // namespace hig
