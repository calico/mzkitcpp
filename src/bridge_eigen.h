#ifndef MZKIT_EIGEN_BRIDGE_H
#define MZKIT_EIGEN_BRIDGE_H

// Force the Eigen core to load and establish the Index type
#include <RcppEigen.h>

// Explicitly pull the Index type into the global or Eigen namespace
// to satisfy the unsupported module solvers
namespace Eigen {
  typedef EIGEN_DEFAULT_DENSE_INDEX_TYPE Index;
}

#endif
