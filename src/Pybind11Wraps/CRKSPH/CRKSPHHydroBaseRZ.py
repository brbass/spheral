#-------------------------------------------------------------------------------
# CRKSPHHydroBaseRZ
#-------------------------------------------------------------------------------
from PYB11Generator import *
from CRKSPHHydroBase import *

@PYB11template()
@PYB11template_dict({"Dimension" : "Dim<2>"})
@PYB11module("SpheralCRKSPH")
class CRKSPHHydroBaseRZ(CRKSPHHydroBase):
    "An area weighted RZ specialization of CRKSPH for cylindrical coordinates"

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    def pyinit(self,
               smoothingScaleMethod = "const SmoothingScaleBase<%(Dimension)s>&",
               Q = "ArtificialViscosity<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               WPi = "const TableKernel<%(Dimension)s>&",
               filter = "const double",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               XSPH = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               correctionOrder = "const RKOrder",
               volumeType = "const RKVolumeType",
               epsTensile = "const double",
               nTensile = "const double",
               limitMultimaterialTopology = "const bool"):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        """Evaluate the derivatives for the principle hydro variables:
mass density, velocity, and specific thermal energy."""
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"
