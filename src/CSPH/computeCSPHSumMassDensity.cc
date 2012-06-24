//------------------------------------------------------------------------------
// Compute the CSPH mass density summation.
//------------------------------------------------------------------------------

#include "computeCSPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
namespace CSPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;

template<typename Dimension>
void
computeCSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& mass,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         const FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // TAU timers.
  TAU_PROFILE("computeCSPHSumMassDensity", "", TAU_USER);

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Zero out the result.
  massDensity = 0.0;

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
    const int firstGhostNodei = nodeList.firstGhostNode();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar mi = mass(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Ai = A(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Self-contribution.
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      massDensity(nodeListi, i) += mi*W0;

      // Get the neighbors for this node (in this NodeList).  We use the approximation here
      // that nodes from other NodeLists do not contribute to the density of this one.
      const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListi, i)[nodeListi];
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;

        // Check if this node pair has already been calculated.
        if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                     nodeListi, j,
                                                     firstGhostNodei)) {
          const Vector& rj = position(nodeListi, j);
          const Scalar& mj = mass(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel weighting and gradient.
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          const Scalar Wi = W.kernelValue(etai, Hdeti);
          const Scalar Wj = W.kernelValue(etaj, Hdetj);

          // Sum the pair-wise contributions.
          massDensity(nodeListi, i) += mj*Wi;
          massDensity(nodeListi, j) += mi*Wj;
        }
      }
      
      // Apply the A correction.
      massDensity(nodeListi, i) *= Ai;
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  namespace CSPHSpace {
    template void computeCSPHSumMassDensity(const ConnectivityMap<Dim<1> >&, 
                                            const TableKernel<Dim<1> >&, 
                                            const FieldList<Dim<1>, Dim<1>::Vector>&,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                            FieldList<Dim<1>, Dim<1>::Scalar>&);
    template void computeCSPHSumMassDensity(const ConnectivityMap<Dim<2> >&, 
                                            const TableKernel<Dim<2> >&, 
                                            const FieldList<Dim<2>, Dim<2>::Vector>&,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                            FieldList<Dim<2>, Dim<2>::Scalar>&);
    template void computeCSPHSumMassDensity(const ConnectivityMap<Dim<3> >&, 
                                            const TableKernel<Dim<3> >&, 
                                            const FieldList<Dim<3>, Dim<3>::Vector>&,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                            FieldList<Dim<3>, Dim<3>::Scalar>&);
  }
}

