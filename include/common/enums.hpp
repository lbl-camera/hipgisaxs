/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: enums.hpp
 *  Created: Jun 11, 2012
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

#ifndef __HIG_ENUMS_HPP__
#define __HIG_ENUMS_HPP__

namespace hig {

  // these are values of string_token:

  enum ShapeName {
    shape_null,
    shape_error,
    shape_box,        /* box */
    shape_cube,       /* same as box */
    shape_cylinder,      /* cylinder */
    shape_sphere,      /* sphere */
    shape_truncpyr,    /* truncated pyramid */
    shape_trunccone,   /* truncated cone */
    shape_prism3,      /* 3-fold prism */
    shape_prism6,      /* 6-fold prism */
    shape_prism3x,     /* triangular grating in x direction */
    shape_sawtooth,    /* sawtooth (prism along x) (same as above) */
    shape_random_cylinders,
    shape_horizontal_cylinder,
    shape_pyramid,
    shape_cone,
    // custom (numerical)
    shape_custom      /* a custom shape */
  }; // enum ShapeName

  enum ShapeParamType {
    param_null,
    param_error,
    param_radius,      /* radius */
    param_xsize,      /* length in x dimension */
    param_ysize,      /* length in y dimension */
    param_height,      /* length in z dimension */
    param_edge,        /* length of an edge (for symmetric objects) */
    param_baseangle      /* base angle */
  }; // enum ShapeParamType

  enum LatticeType {
    lattice_null,
    lattice_error,
    lattice_bcc,      /* bcc */
    lattice_cubic,      /* default cubic */
    lattice_fcc,      /* fcc */
    lattice_fco,      /* fco */
    lattice_hcp,      /* hcp */
    lattice_hex        /* hex */
  }; // enum LatticeType

  // support for liquid crustals
  enum StructureType {
    paracrystal_type,
    percusyevick_type,
    default_type
  };

  enum StatisticType {
    stat_null,
    stat_error,
    stat_none,        /* default */
    stat_cauchy,      /* cauchy/lorentzian distribution */
    stat_gaussian,    /* gaussian distribution */
    stat_random,      /* random distribution */
    stat_range,       /* range of values (basically same as stat_random) */
    stat_t,           /* t distribution */
    stat_uniform      /* uniform distribution */
  }; // enum StatisticType

  enum OutputRegionType {
    region_null,
    region_error,
    region_angles,      /* angles region */
    region_pixels,      /* pixels region */
    region_qspace      /* default qspace region */
  }; // enum OutputRegionType

  enum ShapeFileType {
    shape_file_null,
    shape_file_error,
    shape_file_data,    /* shape file format .dat */
    shape_file_hdf5,    /* shape file in HDF5 format */
    shape_file_object    /* wavefront OBJ object file (e.g. from maya) */
  }; // enum ShapeFileType

  enum StructCorrelationType {
    structcorr_error,    /* error type */
    structcorr_null,    /* default, no correlation */
    structcorr_nGnE,    /* another name for default */
    structcorr_nGE,      /* non-corr grain, corr ensemble */
    structcorr_GnE,      /* corr grain, non-corr ensemble */
    structcorr_GE,      /* both correlated */
  }; // enum StructCorrelationType


  /**
   * fitting related enums
   */

  /**
   * fitting algorithms
   */
  enum FittingAlgorithmName {
    algo_error,         /* error type */
    algo_null,          /* default, no algorithm */
    algo_pounders,      /* pounders algorithm (from tao) */
    algo_pso,           /* particle swarm optimization algorithm */
    algo_lmvm,          /* lmvm algorithm (from tao) */
    algo_bruteforce,    /* brute force optimization */
    algo_none_pounders, /* dont fit, evaluate objective */
    algo_none_lmvm,     /* dont fit, evaluate objective */
    algo_none_pso       /* dont fit, evaluate objective */
  }; // enum FittingAlgorithmName


  /**
   * parameters for various fitting algorithms
   */
  enum FitAlgorithmParamType {
    algo_param_error,           /* error type */
    algo_param_null,            /* default, null parameter */
    algo_pounders_param_delta,  /* delta for pounders algorithm */
    algo_pso_param_omega,       /* omega for pso algorithm */
    algo_pso_param_phi1,        /* phi1 for pso algorithm */
    algo_pso_param_phi2,        /* phi2 for pso algorithm */
    algo_pso_param_nparticle,   /* number of particles for pso algorithm */
    algo_pso_param_ngen,        /* number of generations for pso algorithm */
    algo_pso_param_tune_omega,  /* flag to enable tuning pso omega parameter */
    algo_pso_param_type         /* type of the pso algorithm flavor */
  }; // enum FitAlgorithmParamType


  /**
   * file type of the reference data in fitting
   */
  enum ReferenceFileType {
    reference_file_null,    /* default, null type */
    reference_file_error,    /* error type */
    reference_file_ascii,    /* plain text file */
    reference_file_edf      /* EDF file format */
  }; // ReferenceFileType


  /**
   * distance metrics for fitting reference data comparison
   */
  enum FittingDistanceMetric {
    metric_error,                         /* error type */
    metric_null,                          /* default null type */
    metric_sqrt_unit_norm_l1,             /* sqrt of intensities, unit normalized, l1 distance */
    metric_sqrt_unit_norm_l2,             /* sqrt of intensities, unit normalized, l2 distance */
    metric_sqrt_c_norm_l2,                /* sqrt of intensities, c normalized, l2 distance */
    metric_cbrt_unit_norm_l1,             /* cbrt of intensities, unit normalized, l1 distance */
    metric_cbrt_unit_norm_l2,             /* cbrt of intensities, unit normalized, l2 distance */
    metric_cbrt_c_norm_l2,                /* cbrt of intensities, c normalized, l2 distance */
    metric_sqrt_unit_norm_l1_residual,    /* sqrt of intensities, unit norm, l2 dist, residual */
    metric_sqrt_unit_norm_l2_residual,    /* sqrt of intensities, unit norm, l1 dist, residual */
    metric_sqrt_c_norm_l2_residual        /* sqrt of intensities, c norm, l2 dist, residual */
  }; // enum FittingDistanceMetricType

} // namespace hig

#endif /* _ENUMS_HPP_ */
