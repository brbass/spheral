//---------------------------------Spheral++----------------------------------//
// RKCorrections
//
// Computes RK corrections for other physics packages
//----------------------------------------------------------------------------//
#include "RKCorrections.hh"

#include <limits>
#include "RKUtilities.hh"
#include "computeRKVolumes.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
RKCorrections<Dimension, correctionOrder>::
RKCorrections(const DataBase<Dimension>& dataBase,
              const TableKernel<Dimension>& W,
              const RKVolumeType volumeType,
              const bool needHessian):
  mDataBase(dataBase),
  mW(W),
  mVolumeType(volumeType),
  mNeedHessian(needHessian),
  mVolume(FieldStorageType::CopyFields),
  mZerothCorrections(FieldStorageType::CopyFields),
  mCorrections(FieldStorageType::CopyFields),
  mSurfaceArea(FieldStorageType::CopyFields),
  mNormal(FieldStorageType::CopyFields),
  mSurfacePoint(FieldStorageType::CopyFields),
  mEtaVoidPoints(FieldStorageType::CopyFields),
  mCells(FieldStorageType::CopyFields),
  mCellFaceFlags(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
RKCorrections<Dimension, correctionOrder>::
~RKCorrections() {
}

//------------------------------------------------------------------------------
// Optional hook to initialize once when the problem is starting up
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  // Initialize state
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mCorrections = dataBase.newFluidFieldList(std::vector<double>(), HydroFieldNames::rkCorrections);
  mZerothCorrections = dataBase.newFluidFieldList(std::vector<double>(), HydroFieldNames::rkZerothCorrections);
  mSurfaceArea = dataBase.newFluidFieldList(0.0, HydroFieldNames::surfaceArea);
  mNormal = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::normal);
  
  // Initialize the Voronoi stuff
  mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(std::vector<Vector>(), HydroFieldNames::etaVoidPoints);
  if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
    mCells = dataBase.newFluidFieldList(FacetedVolume(), HydroFieldNames::cells);
    mCellFaceFlags = dataBase.newFluidFieldList(std::vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags);
  }
  mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero, "delta centroid");
  
  // Get some more data
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto mass = dataBase.fluidMass();
  const auto H = dataBase.fluidHfield();
  const auto position = dataBase.fluidPosition();
  const auto massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> damage;
  if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
    damage = dataBase.solidEffectiveDamage();
  }
  
  // Compute the volumes
  computeRKVolumes(connectivityMap, mW,
                   position, mass, massDensity, H, damage,
                   this->boundaryConditions(), mVolumeType,
                   mSurfacePoint, mDeltaCentroid, mEtaVoidPoints, mCells, mCellFaceFlags,
                   mVolume);
  
  // Apply boundaries to newly computed terms
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mVolume);
    if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(mSurfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }
  
  // Compute corrections
  RKUtilities<Dimension, correctionOrder>::
    computeCorrections(connectivityMap, mW, mVolume, position, H,
                       mNeedHessian, mZerothCorrections, mCorrections);

  // Compute normal direction
  RKUtilities<Dimension, RKOrder::ZerothOrder>::
    computeNormal(connectivityMap, mW, mVolume, position, H,
                  mZerothCorrections, mSurfaceArea, mNormal);
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  // Stuff RKCorrections owns
  state.enroll(mVolume);
  state.enroll(mZerothCorrections);
  state.enroll(mCorrections);
  state.enroll(mSurfaceArea);
  state.enroll(mNormal);
  
  state.enroll(mSurfacePoint);
  state.enroll(mEtaVoidPoints);
  if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
    state.enroll(mCells);
    state.enroll(mCellFaceFlags);
  }

  // Stuff RKCorrections needs that might have been enrolled elsewhere
  auto position = dataBase.fluidPosition();
  auto mass = dataBase.fluidMass();
  auto massDensity = dataBase.fluidMassDensity();
  auto H = dataBase.fluidHfield();
  if (not state.registered(position)) state.enroll(position);
  if (not state.registered(mass)) state.enroll(mass);
  if (not state.registered(massDensity)) state.enroll(massDensity);
  if (not state.registered(H)) state.enroll(H);
}

//------------------------------------------------------------------------------
// No derivatives to register
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  // Get state variables
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto zerothCorrections = state.fields(HydroFieldNames::rkZerothCorrections, std::vector<double>());
  auto corrections = state.fields(HydroFieldNames::rkCorrections, std::vector<double>());
  auto surfaceArea = state.fields(HydroFieldNames::surfaceArea, 0.0);
  auto normal = state.fields(HydroFieldNames::normal, Vector::zero);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());

  // Apply ghost boundary conditions
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(vol);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(zerothCorrections);
    (*boundaryItr)->applyFieldListGhostBoundary(corrections);
    (*boundaryItr)->applyFieldListGhostBoundary(surfaceArea);
    (*boundaryItr)->applyFieldListGhostBoundary(normal);
    (*boundaryItr)->applyFieldListGhostBoundary(surfacePoint);
    (*boundaryItr)->applyFieldListGhostBoundary(etaVoidPoints);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // Get state variables
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);

  // Enforce boundary conditions
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(vol);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
  } 
}

//------------------------------------------------------------------------------
// No time step vote
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
typename RKCorrections<Dimension, correctionOrder>::TimeStepType
RKCorrections<Dimension, correctionOrder>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  return std::make_pair(std::numeric_limits<double>::max(), std::string("RKCorrections: no vote"));
}

//------------------------------------------------------------------------------
// Compute new volumes
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // Get data
  const auto& W = mW;
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  FieldList<Dimension, FacetedVolume> cells;
  FieldList<Dimension, std::vector<CellFaceFlag>> cellFaceFlags;
  if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
    cells = state.fields(HydroFieldNames::cells, FacetedVolume());
    cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, std::vector<CellFaceFlag>());
  }
  
  // Compute volumes
  computeRKVolumes(connectivityMap, W,
                   position, mass, massDensity, H, damage,
                   this->boundaryConditions(), mVolumeType,
                   surfacePoint, mDeltaCentroid, mEtaVoidPoints, cells, cellFaceFlags,
                   volume);
  
  // Apply ghost boundaries to Voronoi stuff
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(volume);
    if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(cells);
      (*boundItr)->applyFieldListGhostBoundary(surfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  
}

//------------------------------------------------------------------------------
// Compute new RK corrections
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  // Get data
  const auto& W = mW;
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto zerothCorrections = state.fields(HydroFieldNames::rkZerothCorrections, std::vector<double>());
  auto corrections = state.fields(HydroFieldNames::rkCorrections, std::vector<double>());
  auto surfaceArea = state.fields(HydroFieldNames::surfaceArea, 0.0);
  auto normal = state.fields(HydroFieldNames::normal, Vector::zero);
  
  // Compute corrections
  RKUtilities<Dimension, correctionOrder>::
    computeCorrections(connectivityMap, W, volume, position, H,
                       mNeedHessian, zerothCorrections, corrections);

  RKUtilities<Dimension, RKOrder::ZerothOrder>::
    computeNormal(connectivityMap, W, volume, position, H,
                  zerothCorrections, surfaceArea, normal);
  
  // Apply ghost boundaries to corrections
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(zerothCorrections);
    (*boundaryItr)->applyFieldListGhostBoundary(corrections);
    (*boundaryItr)->applyFieldListGhostBoundary(surfaceArea);
    (*boundaryItr)->applyFieldListGhostBoundary(normal);
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// No derivatives to evaluate
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
}

//------------------------------------------------------------------------------
// Nothing to finalize
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
finalize(const Scalar time, 
         const Scalar dt,
         DataBase<Dimension>& dataBase, 
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Dump the current state to the given file
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
dumpState(FileIO& file, const std::string& pathName) const {
  file.write(mVolume, pathName + "/Volume");
  file.write(mZerothCorrections, pathName + "/ZerothRKCorrections");
  file.write(mCorrections, pathName + "/RKCorrections");
}

//------------------------------------------------------------------------------
// Restore the state from the given file
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
void
RKCorrections<Dimension, correctionOrder>::
restoreState(const FileIO& file, const std::string& pathName) {
  file.read(mVolume, pathName + "/Volume");
  file.read(mZerothCorrections, pathName + "/ZerothRKCorrections");
  file.read(mCorrections, pathName + "/RKCorrections");
}

} // end namespace Spheral
