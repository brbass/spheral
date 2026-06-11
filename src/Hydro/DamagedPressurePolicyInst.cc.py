text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/DamagedPressurePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DamagedPressurePolicy<Dim< %(ndim)s > >;
}
"""
