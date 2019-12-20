text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/gradientRK.cc"
#include "Geometry/Dimension.hh"
#include "SPH/NodeCoupling.hh"

namespace Spheral {
"""

for Value in ("Scalar", "Vector"):
    text += """
template 
FieldList<Dim<%%(ndim)s>, MathTraits<Dim<%%(ndim)s>, Dim<%%(ndim)s>::%(Value)s>::GradientType> 
gradientRK<Dim<%%(ndim)s>, Dim<%%(ndim)s>::%(Value)s>(const FieldList<Dim<%%(ndim)s>, Dim<%%(ndim)s>::%(Value)s>& fieldList,
                                                      const FieldList<Dim<%%(ndim)s>, Dim<%%(ndim)s>::Vector>& position,
                                                      const FieldList<Dim<%%(ndim)s>, Dim<%%(ndim)s>::Scalar>& weight,
                                                      const FieldList<Dim<%%(ndim)s>, Dim<%%(ndim)s>::SymTensor>& H,
                                                      const ConnectivityMap<Dim<%%(ndim)s>>& connectivityMap,
                                                      const TableKernel<Dim<%%(ndim)s>>& W,
                                                      const RKOrder correctionOrder,
                                                      const FieldList<Dim<%%(ndim)s>, std::vector<double>>& corrections,
                                                      const NodeCoupling& nodeCoupling);

""" % {"Value" : Value}

text += """
}
"""