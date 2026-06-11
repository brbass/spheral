//---------------------------------Spheral++----------------------------------//
// DamagedPressurePolicy
// 
// Updates the damaged pressure separately from the pressure to avoid circular
// circular dependencies.
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamagedPressurePolicy_hh__
#define __Spheral_DamagedPressurePolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class DamagedPressurePolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  DamagedPressurePolicy();
  virtual ~DamagedPressurePolicy() = default;

  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  DamagedPressurePolicy(const DamagedPressurePolicy& rhs) = delete;
  DamagedPressurePolicy& operator=(const DamagedPressurePolicy& rhs) = delete;
};

}

#endif
