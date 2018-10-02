#-------------------------------------------------------------------------------
# TensorSVPHViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class TensorSVPHViscosity(ArtificialViscosity):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = ("const Scalar", "1.0"),
               Cquadratic = ("const Scalar", "1.0"),
               fslice = ("const Scalar", "0.5")):
        "TensorSVPHViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    fslice = PYB11property("Scalar", "fslice", "fslice")
    DvDx = PYB11property("const std::vector<Tensor>&", "DvDx", returnpolicy="reference_internal")
    shearCorrection = PYB11property("const std::vector<Scalar>&", "shearCorrection", returnpolicy="reference_internal")
    Qface = PYB11property("const std::vector<Tensor>&", "Qface", returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, TensorSVPHViscosity, virtual=True, pure_virtual=False)
