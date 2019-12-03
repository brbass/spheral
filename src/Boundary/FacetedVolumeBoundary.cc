//---------------------------------Spheral++----------------------------------//
// FacetedVolumeBoundary -- Apply a FacetedVolume boundary condition to Spheral++
// Fields.
//----------------------------------------------------------------------------//
#include "findNodesTouchingThroughPlanes.hh"
#include "mapPositionThroughPlanes.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/innerProduct.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Mesh/Mesh.hh"

#include "FacetedVolumeBoundary.hh"

using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Construct the plane associated with a Facet
//------------------------------------------------------------------------------
// 2D
GeomPlane<Dim<2>>
facetPlane(const GeomFacet2d& facet, const bool interiorBoundary) {
  const auto& p1 = facet.point1();
  const auto  nhat = (interiorBoundary ?
                       facet.normal() :
                      -facet.normal());
  return GeomPlane<Dim<2>>(p1, nhat);
}

//..............................................................................
// 3D
GeomPlane<Dim<3>>
facetPlane(const GeomFacet3d& facet, const bool interiorBoundary) {
  const auto& p1 = facet.point(0);
  const auto  nhat = (interiorBoundary ?
                       facet.normal() :
                      -facet.normal());
  return GeomPlane<Dim<3>>(p1, nhat);
}

//------------------------------------------------------------------------------
// Construct the subvolume associated with a Facet
//------------------------------------------------------------------------------
// 2D
Dim<2>::FacetedVolume
facetSubVolume(const GeomFacet2d& facet, const Dim<2>::Vector& centroid) {
  std::vector<Dim<2>::Vector> points = {facet.point1(), facet.point2(), centroid };
  return Dim<2>::FacetedVolume(points);
}

//..............................................................................
// 3D
Dim<3>::FacetedVolume
facetSubVolume(const GeomFacet3d& facet, const Dim<3>::Vector& centroid) {
  std::vector<Dim<3>::Vector> points = {centroid};
  const auto& ipoints = facet.ipoints();
  const auto  n = ipoints.size();
  for (auto i = 0; i < n; ++i) points.push_back(facet.point(i));
  return Dim<3>::FacetedVolume(points);
}

//------------------------------------------------------------------------------
// Find the node touching a given facet
//------------------------------------------------------------------------------
// Polygon
std::vector<int>
nodesTouchingFacet(const NodeList<Dim<2>>& nodes,
                   const GeomFacet2d& facet,
                   const bool interiorBoundary) {
  const auto& p1 = facet.point1();
  const auto& p2 = facet.point2();
  const auto  phat = (p2 - p1).unitVector();
  const auto  smax = (p2 - p1).magnitude();
  const auto  plane = facetPlane(facet, interiorBoundary);
  const auto potentials = findNodesTouchingThroughPlanes(nodes, plane, plane);
  std::vector<int> result;
  const auto& pos = nodes.positions();
  for (const auto i: potentials) {
    const auto s = (pos(i) - p1).dot(phat);
    if (s >= 0.0 and s <= smax) result.push_back(i);
  }
  return result;
}

//..............................................................................
// Polyhedron
std::vector<int>
nodesTouchingFacet(const NodeList<Dim<3>>& nodes,
                   const GeomFacet3d& facet,
                   const bool interiorBoundary) {
  std::vector<int> result;
  return result;
}

//------------------------------------------------------------------------------
// Copy control->ghost values
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
copyControlValues(Field<Dimension, Value>& field,
                  const std::vector<int>& controls,
                  const std::vector<int>& ghosts) {
  REQUIRE(controls.size() == ghosts.size());
  auto controlItr = controls.begin();
  auto ghostItr = ghosts.begin();
  for (; controlItr < controls.end(); ++controlItr, ++ghostItr) {
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// The per value type method of mapping values
//------------------------------------------------------------------------------
// Vector
template<typename Dimension>
void
mapValue(typename Dimension::Vector& ghost,
         const typename Dimension::Vector& control,
         const typename Dimension::Tensor& R) {
  ghost = R*control;
}

// Tensor
template<typename Dimension>
void
mapValue(typename Dimension::Tensor& ghost,
         const typename Dimension::Tensor& control,
         const typename Dimension::Tensor& R) {
  ghost = R*control*R;
}

// SymTensor
template<typename Dimension>
void
mapValue(typename Dimension::SymTensor& ghost,
         const typename Dimension::SymTensor& control,
         const typename Dimension::Tensor& R) {
  ghost = (R*control*R).Symmetric();
}

// ThirdRankTensor
template<typename Dimension>
void
mapValue(typename Dimension::ThirdRankTensor& ghost,
         const typename Dimension::ThirdRankTensor& control,
         const typename Dimension::Tensor& R) {
  ghost.Zero();
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    for (unsigned j = 0; j != Dimension::nDim; ++j) {
      for (unsigned k = 0; k != Dimension::nDim; ++k) {
        for (unsigned q = 0; q != Dimension::nDim; ++q) {
          for (unsigned r = 0; r != Dimension::nDim; ++r) {
            for (unsigned s = 0; s != Dimension::nDim; ++s) {
              ghost(i,j,k) += R(i,q)*R(j,r)*R(k,s)*control(q,r,s);
            }
          }
        }
      }
    }
  }
}

// FourthRankTensor
template<typename Dimension>
void
mapValue(typename Dimension::FourthRankTensor& ghost,
         const typename Dimension::FourthRankTensor& control,
         const typename Dimension::Tensor& R) {
  ghost.Zero();
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    for (unsigned j = 0; j != Dimension::nDim; ++j) {
      for (unsigned k = 0; k != Dimension::nDim; ++k) {
        for (unsigned l = 0; l != Dimension::nDim; ++l) {
          for (unsigned q = 0; q != Dimension::nDim; ++q) {
            for (unsigned r = 0; r != Dimension::nDim; ++r) {
              for (unsigned s = 0; s != Dimension::nDim; ++s) {
                for (unsigned t = 0; t != Dimension::nDim; ++t) {
                  ghost(i,j,k,l) += R(i,q)*R(j,r)*R(k,s)*R(l,t)*control(q,r,s,t);
                }
              }
            }
          }
        }
      }
    }
  }
}

// FifthRankTensor
template<typename Dimension>
void
mapValue(typename Dimension::FifthRankTensor& ghost,
         const typename Dimension::FifthRankTensor& control,
         const typename Dimension::Tensor& R) {
  ghost.Zero();
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    for (unsigned j = 0; j != Dimension::nDim; ++j) {
      for (unsigned k = 0; k != Dimension::nDim; ++k) {
        for (unsigned l = 0; l != Dimension::nDim; ++l) {
          for (unsigned m = 0; m != Dimension::nDim; ++m) {
            for (unsigned q = 0; q != Dimension::nDim; ++q) {
              for (unsigned r = 0; r != Dimension::nDim; ++r) {
                for (unsigned s = 0; s != Dimension::nDim; ++s) {
                  for (unsigned t = 0; t != Dimension::nDim; ++t) {
                    for (unsigned u = 0; u != Dimension::nDim; ++u) {
                      ghost(i,j,k,l,u) += R(i,q)*R(j,r)*R(k,s)*R(l,t)*R(m,u)*control(q,r,s,t,u);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Map control->ghost values, depending on the Value type
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
mapControlValues(Field<Dimension, Value>& field,
                 const std::vector<typename Dimension::Tensor>& reflectOperators,
                 const std::map<std::string, std::vector<std::vector<int>>>& facetControlNodes,
                 const std::map<std::string, std::vector<std::pair<int,int>>>& facetGhostNodes) {
  const auto& name = field.nodeList().name();
  const auto& controlNodes = facetControlNodes.find(name)->second;
  const auto& ghostNodeRanges = facetGhostNodes.find(name)->second;
  const auto  nfacets = reflectOperators.size();
  REQUIRE(controlNodes.size() == nfacets);
  REQUIRE(ghostNodeRanges.size() == nfacets);
  for (auto ifacet = 0; ifacet < nfacets; ++ifacet) {
    const auto& R = reflectOperators[ifacet];
    const auto& controls = controlNodes[ifacet];
    auto        ghostID = ghostNodeRanges[ifacet].first;
    CHECK(ghostNodeRanges[ifacet].second - ghostID == controls.size());
    for (const auto i: controls) {
      mapValue<Dimension>(field(ghostID), field(i), R);
      ++ghostID;
    }
    CHECK(ghostID == ghostNodeRanges[ifacet].second);
  }
}

//------------------------------------------------------------------------------
// Map control->ghost values, depending on the Value type
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
enforceViolationValues(Field<Dimension, Value>& field,
                       const std::vector<int>& violationNodes,
                       const std::vector<typename Dimension::Tensor>& violationOps) {
  const auto n = violationNodes.size();
  REQUIRE(violationOps.size() == n);
  Value newVal;
  for (auto k = 0; k < n; ++k) {
    const auto i = violationNodes[k];
    mapValue<Dimension>(newVal, field(i), violationOps[k]);
    field(i) = newVal;
  }
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FacetedVolumeBoundary<Dimension>::FacetedVolumeBoundary(const FacetedVolume& poly,
                                                        const bool interiorBoundary,
                                                        const bool useGhosts):
  Boundary<Dimension>(),
  mPoly(poly),
  mInteriorBoundary(interiorBoundary),
  mUseGhosts(useGhosts),
  mReflectOperators(),
  mFacetControlNodes(),
  mFacetGhostNodes(),
  mViolationOperators() {

  // Build the reflection operators for each facet of the poly.
  const auto& facets = poly.facets();
  for (const auto& facet: facets) {
    mReflectOperators.push_back(interiorBoundary ?
                                planarReflectingOperator<Dimension>(facet.normal()) :
                                planarReflectingOperator<Dimension>(-facet.normal()));
  }
  ENSURE(mReflectOperators.size() == facets.size());

}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FacetedVolumeBoundary<Dimension>::~FacetedVolumeBoundary() {
}

//------------------------------------------------------------------------------
// Create ghost nodes for the NodeList
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::setGhostNodes(NodeList<Dimension>& nodeList) {
  if (mUseGhosts) {
    const auto& facets = mPoly.facets();
    const auto  nfacets = facets.size();
    const auto  name = nodeList.name();
    this->addNodeList(nodeList);
    auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
    const auto  centroid = mPoly.centroid();
    auto&       pos = nodeList.positions();
    auto&       H = nodeList.Hfield();

    // Allocate storage for controls and ghosts
    mFacetControlNodes[name] = vector<vector<int>>(nfacets);
    mFacetGhostNodes[name] = vector<pair<int,int>>(nfacets);
    auto& controls = mFacetControlNodes[name];
    auto& ghostRanges = mFacetGhostNodes[name];

    // Build the control points.
    auto firstGhost = nodeList.numNodes();
    vector<Vector> posGhost;
    vector<SymTensor> Hghost;
    for (auto k = 0; k < nfacets; ++k) {
      const auto& facet = facets[k];
      const auto  facetPoly = facetSubVolume(facet, centroid);
      const auto  plane = facetPlane(facets[k], mInteriorBoundary);
      const auto& R = mReflectOperators[k];
      ghostRanges[k].first = firstGhost;
      const auto  potentials = nodesTouchingFacet(nodeList, facet, mInteriorBoundary);
      for (const auto i: potentials) {
        const auto posj = mapPositionThroughPlanes(pos(i), plane, plane);
        if (facetPoly.convexContains(posj, false)) {
          controls[k].push_back(i);
          firstGhost += 1;
          posGhost.push_back(posj);
          Hghost.push_back((R*H(i)*R).Symmetric());
        }
      }
      ghostRanges[k].second = firstGhost;
      CHECK(posGhost.size() == Hghost.size() == (firstGhost - ghostRanges[k].first));
      boundaryNodes.controlNodes.insert(boundaryNodes.controlNodes.end(),
                                        controls[k].begin(),
                                        controls[k].end());
    }
    CHECK(mFacetControlNodes[name].size() == facets.size());
    CHECK(mFacetGhostNodes[name].size() == facets.size());
    CHECK(posGhost.size() == Hghost.size() == boundaryNodes.controlNodes.size())

    // Update the ghost node info.
    const auto numNewGhosts = posGhost.size();
    firstGhost = nodeList.numNodes();
    nodeList.numGhostNodes(nodeList.numGhostNodes() + numNewGhosts);
    std::copy(posGhost.begin(), posGhost.end(), &(pos[firstGhost]));
    std::copy(Hghost.begin(), Hghost.end(), &(H[firstGhost]));
  }
}

//------------------------------------------------------------------------------
// Set the minimal field values necessary to make the ghost nodes valid.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::updateGhostNodes(NodeList<Dimension>& nodeList) {

  const auto  name = nodeList.name();
  const auto& controls = mFacetControlNodes[name];
  const auto& ghostRanges = mFacetGhostNodes[name];
  const auto& facets = mPoly.facets();
  const auto  nfacets = facets.size();
  auto&       pos = nodeList.positions();
  auto&       H = nodeList.Hfield();
  CHECK(controls.size() == nfacets);
  CHECK(ghostRanges.size() == nfacets);

  // Walk each of the facets
  for (auto ifacet = 0; ifacet < nfacets; ++ifacet) {
    const auto  n = controls[ifacet].size();
    const auto  plane = facetPlane(facets[ifacet], mInteriorBoundary);
    const auto& R = mReflectOperators[ifacet];
    CHECK(ghostRanges[ifacet].second - ghostRanges[ifacet].first == n);
    for (auto k = 0; k < n; ++k) {
      const auto i = controls[ifacet][k];
      const auto j = ghostRanges[ifacet].first + k;
      pos(j) = mapPositionThroughPlanes(pos(i), plane, plane);
      H(j) = (R*H(i)*R).Symmetric();
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, just perform a copy.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
  const auto& nodeList = field.nodeList();
  copyControlValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for scalar fields, just perform a copy.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  const auto& nodeList = field.nodeList();
  copyControlValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for Vector fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);}

// Specialization for Tensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
}

// Specialization for ThirdRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
}

// Specialization for FourthRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
}

// Specialization for FifthRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
}

// Specialization for FacetedVolumes
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  const auto& nodeList = field.nodeList();
  copyControlValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));   // Punt for FacetedVolumes and just copy them for now
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::setViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  this->addNodeList(nodeList);
  auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& vNodes = boundaryNodes.violationNodes;
  vNodes.clear();

  // Loop over all the internal nodes in the NodeList and look for any that have
  // wandered into the excluded region.
  const auto& pos = nodeList.positions();
  const auto  n = nodeList.numInternalNodes();
  for (auto i = 0; i < n; ++i) {
    const auto interior = mPoly.contains(pos(i), false);
    if (mInteriorBoundary) {
      if (interior) vNodes.push_back(i);
    } else {
      if (not interior) vNodes.push_back(i);
    }
  }

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::updateViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the set of violation nodes for this NodeList.
  const auto& vNodes = this->violationNodes(nodeList);
  const auto  numViolation = vNodes.size();
  auto&       violationOps = mViolationOperators[nodeList.name()];
  violationOps = vector<Tensor>(vNodes.size(), Tensor::one);
  const auto& facets = mPoly.facets();

  // Find the longest scale in the FacetedVolume
  const auto chordLength = (mPoly.xmax() - mPoly.xmin()).magnitude();

  // Loop over these nodes, and reset their positions to valid values.
  auto&       pos = nodeList.positions();
  auto&       H = nodeList.Hfield();
  const auto& vel = nodeList.velocity();
  vector<unsigned> potentialFacets;
  vector<Vector> intersections;
  Vector newPos, newVel;
  SymTensor newH;
  bool inViolation;
  for (auto k = 0; k < numViolation; ++k) {
    auto  i = vNodes[k];
    auto& R = violationOps[k];
    newPos = pos(i);
    newVel = vel(i);
    while (inViolation) {
      // Backtrack to which facet we think the point passed through.
      const auto backPos = newPos - chordLength*newVel.unitVector();
      mPoly.intersect(backPos, newPos, potentialFacets, intersections);
      CHECK(potentialFacets.size() == intersections.size());
      auto minFacet = potentialFacets[0];
      auto minDist = (intersections[0] - newPos).magnitude2();
      for (auto kk = 1; kk < intersections.size(); ++kk) {
        if ((intersections[kk] - newPos).magnitude2() < minDist) {
          minFacet = potentialFacets[kk];
          minDist = (intersections[kk] - newPos).magnitude2();
        }
      }
      const auto plane = facetPlane(facets[minFacet], mInteriorBoundary);
      R *= mReflectOperators[minFacet];
      newPos = mapPositionThroughPlanes(newPos, plane, plane);
      newVel = mReflectOperators[minFacet]*newVel;
      inViolation = mPoly.contains(newPos, false);
      if (mInteriorBoundary) inViolation = not inViolation;
    }
    pos(i) = newPos;
  }

  // Set the Hfield.
  auto& Hfield = nodeList.Hfield();
  this->enforceBoundary(Hfield);
}    

//------------------------------------------------------------------------------
// Enforce the boundary condition on the set of nodes in violation of the 
// boundary.
//------------------------------------------------------------------------------
// Specialization for int fields.  A no-op.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields.  A no-op.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for vector fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for third rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for fourth rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for fifth rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for FacetedVolumes
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  // Punt on FacetedVolumes for now.
}

//------------------------------------------------------------------------------
// Clear out any NodeList information that is currently present.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::reset(const DataBase<Dimension>& dataBase) {
  Boundary<Dimension>::reset(dataBase);
  mFacetControlNodes.clear();
  mFacetGhostNodes.clear();
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#ifdef SPHERAL2D
template class FacetedVolumeBoundary<Dim<2>>;
#endif

#ifdef SPHERAL3D
template class FacetedVolumeBoundary<Dim<3>>;
#endif

}