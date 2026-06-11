//---------------------------------Spheral++----------------------------------//
// DamagedPressurePolicy
// 
// Updates the damaged pressure separately from the pressure to avoid circular
// circular dependencies.
//----------------------------------------------------------------------------//

#include "DamagedPressurePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "Utilities/DBC.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamagedPressurePolicy<Dimension>::
DamagedPressurePolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::pressure,
                                        SolidFieldNames::porosityAlpha,
                                        SolidFieldNames::tensorDamage}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamagedPressurePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::damagedPressure);
  auto& Pd = state.field(key, Scalar());

  // The (undamaged, bulk) pressure this damaged pressure is derived from.
  const auto buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);

  // No damage for this material -> damaged pressure is just the pressure.
  if (not state.registered(buildKey(SolidFieldNames::tensorDamage))) {
    Pd = P;
    return;
  }
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero());

  // Grab the eos for the minimum-pressure thresholds
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(Pd.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();
  const auto* solidEOSptr = dynamic_cast<const SolidEquationOfState<Dimension>*>(&eos);
  const auto Pmin = eos.minimumPressure();
  const auto PminDamage = (solidEOSptr != nullptr ?
                           solidEOSptr->minimumPressureDamage() :
                           Pmin);

  // With porosity the pressure field holds the bulk value (P_solid/alpha).  The
  // damage scaling is defined on the solid pressure, so reconstruct it as
  // P_solid = P*alpha, scale by damage, then divide back by alpha.  Without
  // porosity alpha == 1 and this reduces to operating directly on P.
  const auto usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));
  const Field<Dimension, Scalar>* alphaPtr =
    (usePorosity ? &state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0) : nullptr);

  const auto ni = Pd.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {
    const auto alphai = (usePorosity ? (*alphaPtr)(i) : 1.0);
    const auto Psolid = P(i) * alphai;
    const auto Di = std::max(0.0, std::min(1.0, D(i).eigenValues().maxElement()));
    CHECK(Di >= 0.0 and Di <= 1.0);
    const auto Pdsolid = std::max(Pmin, (1.0 - Di)*Psolid) + std::max(PminDamage, Di*Psolid);
    Pd(i) = Pdsolid / alphai;
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DamagedPressurePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  const auto rhsPtr = dynamic_cast<const DamagedPressurePolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}
