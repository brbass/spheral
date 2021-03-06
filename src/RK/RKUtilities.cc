//---------------------------------Spheral++----------------------------------//
// RKUtilities
//
// Computes and evaluates RK corrections
//----------------------------------------------------------------------------//

#include "RKUtilities.hh"
#include "Eigen/Dense"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Evaluate the base kernel value, gradient, or hessian
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
typename Dimension::Scalar
RKUtilities<Dimension, correctionOrder>::
evaluateBaseKernel(const TableKernel<Dimension>& kernel,
                   const Vector& x,
                   const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto Hdet = H.Determinant();
  return kernel.kernelValue(etaMag, Hdet);
}

template<typename Dimension, RKOrder correctionOrder>
typename Dimension::Vector
RKUtilities<Dimension, correctionOrder>::
evaluateBaseGradient(const TableKernel<Dimension>& kernel,
                     const Vector& x,
                     const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto etaUnit = eta.unitVector();
  const auto Hdet = H.Determinant();
  const auto dk = kernel.gradValue(etaMag, Hdet);
  const auto HetaUnit = H * etaUnit;
  return HetaUnit * dk;
}

template<typename Dimension, RKOrder correctionOrder>
typename Dimension::SymTensor
RKUtilities<Dimension, correctionOrder>::
evaluateBaseHessian(const TableKernel<Dimension>& kernel,
                    const Vector& x,
                    const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto etaUnit = eta.unitVector();
  const auto etaMagInv = safeInv(etaMag);
  const auto Hdet = H.Determinant();
  const auto dk = kernel.gradValue(etaMag, Hdet);
  const auto ddk = kernel.grad2Value(etaMag, Hdet);
  const auto HetaUnit = H * etaUnit;
  const auto HetaUnit2 = HetaUnit.selfdyad();
  const auto H2 = H.square();
  return (H2 - HetaUnit2) * etaMagInv * dk + HetaUnit2 * ddk;
}

template<typename Dimension, RKOrder correctionOrder>
std::pair<typename Dimension::Scalar, typename Dimension::Vector>
RKUtilities<Dimension, correctionOrder>::
evaluateBaseKernelAndGradient(const TableKernel<Dimension>& kernel,
                              const Vector& x,
                              const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto etaUnit = eta.unitVector();
  const auto Hdet = H.Determinant();
  const auto k = kernel.kernelValue(etaMag, Hdet);
  const auto dk = kernel.gradValue(etaMag, Hdet);
  const auto HetaUnit = H * etaUnit;
  return std::make_pair(k, HetaUnit * dk);
}

//------------------------------------------------------------------------------
// Evaluate the RK kernel, gradient, or Hessian
// W^{R}=C^{\top}PW
// W_{\gamma}^{R}=\left(C^{\top}P_{\gamma}+C_{\gamma}^{\top}P\right)W+C^{\top}PW_{\gamma}
// W_{\gamma\zeta}^{R}=\left(C^{\top}P_{\gamma\zeta}+C_{\gamma}^{\top}P_{\zeta}+C_{\zeta}^{\top}P_{\gamma}+C_{\gamma\zeta}^{\top}P\right)W+\left(C^{\top}P_{\gamma}+C_{\gamma}^{\top}P\right)W_{\zeta}+\left(C^{\top}P_{\zeta}+C_{\zeta}^{\top}P\right)W_{\gamma}+C^{\top}PW_{\gamma\zeta}
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
typename Dimension::Scalar
RKUtilities<Dimension, correctionOrder>::
evaluateKernel(const TableKernel<Dimension>& kernel,
               const Vector& x,
               const SymTensor& H,
               const std::vector<double>& corrections) {
  CHECK(corrections.size() == correctionsSize(false)
        || corrections.size() == correctionsSize(true));
  
  // Get kernel and polynomials
  const auto w = evaluateBaseKernel(kernel, x, H);
  PolyArray P;
  getPolynomials(x, P);
  
  // Perform inner products and return result
  const auto CP = innerProductRK(corrections, P, 0, 0);
  return CP * w;
}

template<typename Dimension, RKOrder correctionOrder>
typename Dimension::Vector
RKUtilities<Dimension, correctionOrder>::
evaluateGradient(const TableKernel<Dimension>& kernel,
                 const Vector& x,
                 const SymTensor& H,
                 const std::vector<double>& corrections) {
  CHECK(corrections.size() == correctionsSize(false)
        || corrections.size() == correctionsSize(true));
  
  // Get kernel and polynomials
  const auto wdw = evaluateBaseKernelAndGradient(kernel, x, H);
  const auto w = wdw.first;
  const auto dw = wdw.second;
  PolyArray P;
  GradPolyArray dP;
  getPolynomials(x, P);
  getGradPolynomials(x, dP);
  
  // Perform inner products and return result
  const auto CP = innerProductRK(corrections, P, 0, 0);
  Vector result = Vector::zero;
  for (auto d = 0; d < Dimension::nDim; ++ d) {
    const auto CdP = innerProductRK(corrections, dP, 0, offsetGradP(d));
    const auto dCP = innerProductRK(corrections, P, offsetGradC(d), 0);
    result(d) = (CdP + dCP) * w + CP * dw(d);
  }
  return result;
}

template<typename Dimension, RKOrder correctionOrder>
typename Dimension::SymTensor
RKUtilities<Dimension, correctionOrder>::
evaluateHessian(const TableKernel<Dimension>& kernel,
                const Vector& x,
                const SymTensor& H,
                const std::vector<double>& corrections) {
  CHECK(corrections.size() == correctionsSize(false)
        || corrections.size() == correctionsSize(true));
  
  // Get kernel and polynomials
  const auto wdw = evaluateBaseKernelAndGradient(kernel, x, H);
  const auto w = wdw.first;
  const auto dw = wdw.second;
  const auto ddw = evaluateBaseHessian(kernel, x, H);
  PolyArray P;
  GradPolyArray dP;
  HessPolyArray ddP;
  getPolynomials(x, P);
  getGradPolynomials(x, dP);
  getHessPolynomials(x, ddP);
  
  // Perform inner products and return result
  // Could precompute the inner products in a separate loop for efficiency
  const auto CP = innerProductRK(corrections, P, 0, 0);
  SymTensor result = SymTensor::zero;
  for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
    const auto Cd1P = innerProductRK(corrections, dP, 0, offsetGradP(d1));
    const auto d1CP = innerProductRK(corrections, P, offsetGradC(d1), 0);
    for (auto d2 = d1; d2 < Dimension::nDim; ++d2) {
      const auto Cd2P = innerProductRK(corrections, dP, 0, offsetGradP(d2));
      const auto d2CP = innerProductRK(corrections, P, offsetGradC(d2), 0);
      const auto CddP = innerProductRK(corrections, ddP, 0, offsetHessP(d1, d2));
      const auto d1Cd2P = innerProductRK(corrections, dP, offsetGradC(d1), offsetGradP(d2));
      const auto d2Cd1P = innerProductRK(corrections, dP, offsetGradC(d2), offsetGradP(d1));
      const auto ddCP = innerProductRK(corrections, P, offsetHessC(d1, d2), 0);
      
      result(d1, d2) = (CddP + d1Cd2P + d2Cd1P + ddCP) * w + (Cd1P + d1CP) * dw(d2) + (Cd2P + d2CP) * dw(d1) + CP * ddw(d1, d2);
    }
  }
  return result;
}

template<typename Dimension, RKOrder correctionOrder>
std::pair<typename Dimension::Scalar, typename Dimension::Vector>
RKUtilities<Dimension, correctionOrder>::
evaluateKernelAndGradient(const TableKernel<Dimension>& kernel,
                          const Vector& x,
                          const SymTensor& H,
                          const std::vector<double>& corrections) {
  CHECK(corrections.size() == correctionsSize(false)
        || corrections.size() == correctionsSize(true));
  
  // Get kernel and polynomials
  const auto wdw = evaluateBaseKernelAndGradient(kernel, x, H);
  const auto w = wdw.first;
  const auto dw = wdw.second;
  PolyArray P;
  GradPolyArray dP;
  getPolynomials(x, P);
  getGradPolynomials(x, dP);
  
  // Perform inner products and return result
  const auto CP = innerProductRK(corrections, P, 0, 0);
  Vector result = Vector::zero;
  for (auto d = 0; d < Dimension::nDim; ++ d) {
    const auto CdP = innerProductRK(corrections, dP, 0, offsetGradP(d));
    const auto dCP = innerProductRK(corrections, P, offsetGradC(d), 0);
    result(d) = (CdP + dCP) * w + CP * dw(d);
  }
  return std::make_pair(CP * w, result);
}

//------------------------------------------------------------------------------
// Compute the corrections
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKUtilities<Dimension, correctionOrder>::
computeCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                   const TableKernel<Dimension>& kernel,
                   const FieldList<Dimension, Scalar>& volume,
                   const FieldList<Dimension, Vector>& position,
                   const FieldList<Dimension, SymTensor>& H,
                   const bool needHessian,
                   FieldList<Dimension, std::vector<double>>& zerothCorrections, 
                   FieldList<Dimension, std::vector<double>>& corrections) {
  // Typedefs: Eigen requires aligned allocator for stl containers before c++17
  typedef Eigen::Matrix<double, polynomialSize, 1> VectorType;
  typedef Eigen::Matrix<double, polynomialSize, polynomialSize> MatrixType;
  typedef std::vector<VectorType, Eigen::aligned_allocator<VectorType>> VectorOfVectorType;
  typedef std::vector<MatrixType, Eigen::aligned_allocator<MatrixType>> VectorOfMatrixType;
  
  // Size info
  const auto numNodeLists = volume.size();
  const auto hessSize = needHessian ? hessBaseSize : 0;
  const auto zerothCorrSize = zerothCorrectionsSize(needHessian);
  const auto corrSize = correctionsSize(needHessian);

  // Check things
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(corrections.size() == numNodeLists);
  
  // Initialize Eigen stuff
  MatrixType M;
  VectorOfMatrixType dM(Dimension::nDim);
  VectorOfMatrixType ddM(hessSize);
  VectorType C;
  VectorOfVectorType dC(Dimension::nDim);
  VectorOfVectorType ddC(hessSize);
  VectorType rhs;

  // Initialize polynomial arrays
  PolyArray p;
  GradPolyArray dp;
  HessPolyArray ddp;

  // Get function for adding contribution to matrices
  auto addToMatrix = [&](const int nodeListi,
                         const int nodei,
                         const int nodeListj,
                         const int nodej) {
    // Get data for point i
    const auto xi = position(nodeListi , nodei);
    
    // Get data for point j
    const auto xj = position(nodeListj, nodej);
    const auto xij = xi - xj;
    const auto Hj = H(nodeListj, nodej);
    const auto vj = volume(nodeListj, nodej);
    const auto wdw = evaluateBaseKernelAndGradient(kernel, xij, Hj);
    
    // Add to matrix
    const auto w = wdw.first;
    getPolynomials(xij, p);
    for (auto k = 0; k < polynomialSize; ++k) {
      for (auto l = k; l < polynomialSize; ++l) {
        M(k, l) += vj * p[k] * p[l] * w;
      }
    }
    
    // Add to gradient matrix
    const auto dw = wdw.second;
    getGradPolynomials(xij, dp);
    for (auto d = 0; d < Dimension::nDim; ++d) {
      const auto offd = offsetGradP(d);
      for (auto k = 0; k < polynomialSize; ++k) {
        for (auto l = k; l < polynomialSize; ++l) {
          dM[d](k,l) += vj * ((dp[offd+k] * p[l] + p[k] * dp[offd+l]) * w + p[k] * p[l] * dw(d));
        }
      }
    }
    
    // Add to Hessian matrix
    if (needHessian) {
      const auto ddw = evaluateBaseHessian(kernel, xij, Hj);
      getHessPolynomials(xij, ddp);
      for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
        const auto offd1 = offsetGradP(d1);
        for (auto d2 = d1; d2 < Dimension::nDim; ++d2) {
          const auto offd2 = offsetGradP(d2);
          const auto offd12 = offsetHessP(d1, d2);
          const auto d12 = flatSymmetricIndex(d1, d2);
          for (auto k = 0; k < polynomialSize; ++k) {
            for (auto l = k; l < polynomialSize; ++l) {
              ddM[d12](k,l) += vj * ((ddp[offd12+k] * p[l] + dp[offd1+k] * dp[offd2+l] + dp[offd2+k] * dp[offd1+l] + p[k] * ddp[offd12+l]) * w + (dp[offd1+k] * p[l] + p[k] * dp[offd1+l]) * dw(d2) + (dp[offd2+k] * p[l] + p[k] * dp[offd2+l]) * dw(d1) + p[k] * p[l] * ddw(d1, d2));
            }
          }
        }
      }
    }
    
    return;
  };
  
  // Compute corrections for each point
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = connectivityMap.numNodes(nodeListi);
    for (auto nodei = 0; nodei < numNodes; ++nodei) {
      // Initialize polynomial matrices for point i
      M.setZero();
      for (auto& mat : dM) {
        mat.setZero();
      }
      for (auto& mat : ddM) {
        mat.setZero();
      }
                            
      // Add contribution from other points
      const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, nodei);
      for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
        for (auto nodej : connectivity[nodeListj]) {
          addToMatrix(nodeListi, nodei, nodeListj, nodej);
        } // nodej
      } // nodeListj

      // Add self contribution
      addToMatrix(nodeListi, nodei, nodeListi, nodei);
      
      // M symmetries
      for (auto k = 0; k < polynomialSize; ++k) {
        for (auto l = 0; l < k; ++l) {
          M(k, l) = M(l, k);
        }
      }

      // dM symmetries
      for (auto d = 0; d < Dimension::nDim; ++d) {
        for (auto k = 0; k < polynomialSize; ++k) {
          for (auto l = 0; l < k; ++l) {
            dM[d](k, l) = dM[d](l, k);
          }
        }
      }

      // ddM symmetries
      if (needHessian) {
        for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
          for (auto d2 = d1; d2 < Dimension::nDim; ++d2) {
            const auto d12 = flatSymmetricIndex(d1, d2);
            for (auto k = 0; k < polynomialSize; ++k) {
              for (auto l = 0; l < k; ++l) {
                ddM[d12](k,l) = ddM[d12](l,k);
              }
            }
          }
        }
      }
      
      // Get inverse of M matrix
      auto solver = M.colPivHouseholderQr();
      
      // Compute corrections
      rhs(0) = 1;
      for (auto k = 1; k < polynomialSize; ++k) {
        rhs(k) = 0;
      }
      C = solver.solve(rhs);

      // Compute gradient corrections
      for (auto d = 0; d < Dimension::nDim; ++d) {
        rhs = -(dM[d] * C);
        dC[d] = solver.solve(rhs);
      }
      
      // Compute hessian corrections
      if (needHessian) {
        for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
          for (auto d2 = d1; d2 < Dimension::nDim; ++d2) {
            const auto d12 = flatSymmetricIndex(d1, d2);
            rhs = -(ddM[d12] * C + dM[d1] * dC[d2] + dM[d2] * dC[d1]);
            ddC[d12] = solver.solve(rhs);
          }
        }
      }

      // Initialize corrections vector
      auto& corr = corrections(nodeListi, nodei);
      corr.resize(corrSize);
      
      // Put corrections into vector
      for (auto k = 0; k < polynomialSize; ++k) {
        corr[k] = C(k);
      }

      // Put gradient corrections into vector
      for (auto d = 0; d < Dimension::nDim; ++d) {
        const auto offd = offsetGradC(d);
        for (auto k = 0; k < polynomialSize; ++k) {
          corr[offd+k] = dC[d](k);
        }
      }

      // Put hessian corrections into vector
      if (needHessian) {
        for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
          for (auto d2 = d1; d2 < Dimension::nDim; ++d2) {
            const auto d12 = flatSymmetricIndex(d1, d2);
            const auto offd12 = offsetHessC(d1, d2);
            for (auto k = 0; k < polynomialSize; ++k) {
              corr[offd12+k] = ddC[d12](k);
            }
          }
        }
      }

      // Initialize zeroth corrections vector
      auto& zerothCorr = zerothCorrections(nodeListi, nodei);
      zerothCorr.resize(zerothCorrSize);
      
      // Compute zeroth correction
      const auto C0 = safeInv(M(0,0));
      zerothCorr[0] = C0;

      // Compute zeroth gradient
      for (auto d = 0; d < Dimension::nDim; ++d) {
        const auto offd = RKUtilities<Dimension, RKOrder::ZerothOrder>::offsetGradC(d);
        zerothCorr[offd] = -dM[d](0,0) * C0 * C0;
      }
      
      // Compute zeroth hessian
      if (needHessian) {
        for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
          const auto offd1 = RKUtilities<Dimension, RKOrder::ZerothOrder>::offsetGradC(d1);
          const auto C1 = zerothCorr[offd1];
          for (auto d2 = d1; d2 < Dimension::nDim; ++d2) {
            const auto d12 = flatSymmetricIndex(d1, d2);
            const auto offd2 = RKUtilities<Dimension, RKOrder::ZerothOrder>::offsetGradC(d2);
            const auto offd12 = RKUtilities<Dimension, RKOrder::ZerothOrder>::offsetHessC(d1, d2);
            const auto C2 = zerothCorr[offd2];
            zerothCorr[offd12] = -(ddM[d12](0,0) * C0 + dM[d1](0,0) * C2 + dM[d2](0,0) * C1) * C0;
          }
        }
      }
    } // nodei
  } // nodeListi
} // computeCorrections

//------------------------------------------------------------------------------
// Compute a guess for the normal direction
// S_{i}n_{i}^{\alpha}=V_{i}\dfrac{\sum_{j}V_{j}\left[\partial_{x_{j}}^{\alpha}\psi_{i}\left(x_{j}\right)+\partial_{x_{i}}^{\alpha}\psi_{j}\left(x_{i}\right)\right]}{\sum_{j}V_{j}\psi_{j}\left(x_{i}\right)}
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKUtilities<Dimension, correctionOrder>::
computeNormal(const ConnectivityMap<Dimension>& connectivityMap,
              const TableKernel<Dimension>& kernel,
              const FieldList<Dimension, Scalar>& volume,
              const FieldList<Dimension, Vector>& position,
              const FieldList<Dimension, SymTensor>& H,
              const FieldList<Dimension, std::vector<double>>& corrections,
              FieldList<Dimension, Scalar>& surfaceArea,
              FieldList<Dimension, Vector>& normal) {
  // Size info
  const auto numNodeLists = volume.size();
  
  // Check things
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(corrections.size() == numNodeLists);
  REQUIRE(surfaceArea.size() == numNodeLists);
  REQUIRE(normal.size() == numNodeLists);
  
  // Get function for adding contribution to surface area and normal
  auto addToNormal = [&](const int nodeListi,
                         const int nodei,
                         const int nodeListj,
                         const int nodej) {
    // Get data for point i
    const auto xi = position(nodeListi , nodei);
    const auto Hi = H(nodeListi, nodei);
    const auto& Ci = corrections(nodeListi, nodei);

    // Get data for point j
    const auto xj = position(nodeListj, nodej);
    const auto Hj = H(nodeListj, nodej);
    const auto vj = volume(nodeListj, nodej);
    const auto& Cj = corrections(nodeListj, nodej);

    // Get kernels
    const auto xij = xi - xj;
    const auto xji = xj - xi;
    const auto wdwij = evaluateKernelAndGradient(kernel, xij, Hj, Ci);
    const auto wij = wdwij.first;
    const auto dwij = wdwij.second;
    const auto dwji = evaluateGradient(kernel, xji, Hi, Cj);
    
    // Add denominator value to surface area
    surfaceArea(nodeListi, nodei) += vj * wij;
    
    // Add numerator value to normal
    normal(nodeListi, nodei) += vj * (dwij + dwji);
  };

  // Sum up normal directions
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = connectivityMap.numNodes(nodeListi);
    for (auto nodei = 0; nodei < numNodes; ++nodei) {
      // Zero out normal
      normal(nodeListi, nodei) = Vector::zero;
      
      // Add contribution from other points
      const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, nodei);
      for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
        for (auto nodej : connectivity[nodeListj]) {
          addToNormal(nodeListi, nodei, nodeListj, nodej);
        }
      }
                            
      // Add self contribution
      addToNormal(nodeListi, nodei, nodeListi, nodei);

      // Divide the numerator by the denominator to get S_i n_i^\alpha
      normal(nodeListi, nodei) *= volume(nodeListi, nodei) / surfaceArea(nodeListi, nodei);
      
      // Get the surface area
      surfaceArea(nodeListi, nodei) = normal(nodeListi, nodei).magnitude();
      
      // Get the normal direction
      normal(nodeListi, nodei) = normal(nodeListi, nodei).unitVector();
    }
  }
}
  
} // end namespace Spheral
