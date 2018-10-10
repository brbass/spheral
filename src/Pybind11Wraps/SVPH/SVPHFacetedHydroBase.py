#-------------------------------------------------------------------------------
# SVPHFactedHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *

from GenericHydro import *
from RestartMethods import *

@PYB11template("Dimension")
class SVPHFacetedHydroBase(GenericHydro):
    "SVPHFacetedHydroBase -- The fluid SVPH faceted hydro algorithm"

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef typename DIM::FourthRankTensor FourthRankTensor;
    typedef typename DIM::FifthRankTensor FifthRankTensor;
    typedef typename Physics<DIM>::TimeStepType TimeStepType;
"""

    def pyinit(self,
               smoothingScaleMethod = "const SmoothingScaleBase<DIM>&",
               W = "const TableKernel<DIM>&",
               Q = "ArtificialViscosity<DIM>&",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               XSVPH = "const bool",
               linearConsistent = "const bool",
               generateVoid = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               fcentroidal = "const Scalar",
               fcellPressure = "const Scalar",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SVPH Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, dataBase = "DataBase<DIM>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(self,
           dataBase = "const DataBase<DIM>&",
           state = "const State<DIM>&",
           derivs = "const StateDerivatives<DIM>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<DIM>&",
                      state = "State<DIM>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<DIM>&",
                   state = "State<DIM>&",
                   derivs = "StateDerivatives<DIM>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "void"
                          
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        """Evaluate the derivatives for the principle hydro variables:
mass density, velocity, and specific thermal energy."""
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<DIM>&",
                 state = "State<DIM>&",
                 derivs = "StateDerivatives<DIM>&"):
        "Finalize the hydro at the completion of an integration step."
        return "void"
                  
    @PYB11virtual
    @PYB11const
    def requireConnectivity(self):
        "This algorithm does not use node->node connectivity."
        return "bool"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<DIM>&",
                             derivs = "StateDerivatives<DIM>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    #...........................................................................
    # Properties
    densityUpdate = PYB11property("MassDensityType", "densityUpdate", "densityUpdate",
                                  doc="Flag to choose whether we want to sum for density, or integrate the continuity equation.")
    HEvolution = PYB11property("HEvolutionType", "HEvolution", "HEvolution",
                               doc="Flag to select how we want to evolve the H tensor.")
    compatibleEnergyEvolution = PYB11property("bool", "compatibleEnergyEvolution", "compatibleEnergyEvolution",
                                              doc="Flag to determine if we're using the total energy conserving compatible energy evolution scheme.")
    XSVPH = PYB11property("bool", "XSVPH", "XSVPH",
                          doc="Flag to determine if we're using the XSVPH algorithm.")
    linearConsistent = PYB11property("bool", "linearConsistent", "linearConsistent",
                                     doc="Flag to select whether or not to use the linear corrections.")
    generateVoid = PYB11property("bool", "generateVoid", "generateVoid",
                                 doc="Flag to select whether or not to generate void points in the tessellation.")
    fcentroidal = PYB11property("Scalar", "fcentroidal", "fcentroidal",
                                doc="Fraction of centroidal motion to apply each step.")
    fcellPressure = PYB11property("Scalar", "fcellPressure", "fcellPressure",
                                  doc="Fraction of the pressure to take from local cell.")
    xmin = PYB11property("const Vector&", "xmin", "xmin",
                         doc="Optionally we can provide a bounding box for use generating the mesh.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",
                         doc="Optionally we can provide a bounding box for use generating the mesh.")
    smoothingScaleMethod = PYB11property("const SmoothingScaleBase<DIM>&", "smoothingScaleMethod",
                                         doc="The object defining how we evolve smoothing scales.")
    mesh = PYB11property("Mesh<DIM>&", "mesh",
                         doc="The tessellation")


    timeStepMask = PYB11property("const FieldList<DIM, int>&", "timeStepMask", returnpolicy="reference_internal")
    pressure = PYB11property("const FieldList<DIM, Scalar>&", "pressure", returnpolicy="reference_internal")
    cellPressure = PYB11property("const FieldList<DIM, Scalar>&", "cellPressure", returnpolicy="reference_internal")
    soundSpeed = PYB11property("const FieldList<DIM, Scalar>&", "soundSpeed", returnpolicy="reference_internal")
    volume = PYB11property("const FieldList<DIM, Scalar>&", "volume", returnpolicy="reference_internal")
    specificThermalEnergy0 = PYB11property("const FieldList<DIM, Scalar>&", "specificThermalEnergy0", returnpolicy="reference_internal")
    Hideal = PYB11property("const FieldList<DIM, SymTensor>&", "Hideal", returnpolicy="reference_internal")
    maxViscousPressure = PYB11property("const FieldList<DIM, Scalar>&", "maxViscousPressure", returnpolicy="reference_internal")
    massDensitySum = PYB11property("const FieldList<DIM, Scalar>&", "massDensitySum", returnpolicy="reference_internal")
    weightedNeighborSum = PYB11property("const FieldList<DIM, Scalar>&", "weightedNeighborSum", returnpolicy="reference_internal")
    massSecondMoment = PYB11property("const FieldList<DIM, SymTensor>&", "massSecondMoment", returnpolicy="reference_internal")
    XSVPHDeltaV = PYB11property("const FieldList<DIM, Vector>&", "XSVPHDeltaV", returnpolicy="reference_internal")
    DxDt = PYB11property("const FieldList<DIM, Vector>&", "DxDt", returnpolicy="reference_internal")
    DvDt = PYB11property("const FieldList<DIM, Vector>&", "DvDt", returnpolicy="reference_internal")
    DmassDensityDt = PYB11property("const FieldList<DIM, Scalar>&", "DmassDensityDt", returnpolicy="reference_internal")
    DspecificThermalEnergyDt = PYB11property("const FieldList<DIM, Scalar>&", "DspecificThermalEnergyDt", returnpolicy="reference_internal")
    DHDt = PYB11property("const FieldList<DIM, SymTensor>&", "DHDt", returnpolicy="reference_internal")
    DvDx = PYB11property("const FieldList<DIM, Tensor>&", "DvDx", returnpolicy="reference_internal")
    internalDvDx = PYB11property("const FieldList<DIM, Tensor>&", "internalDvDx", returnpolicy="reference_internal")
    faceForce = PYB11property("const FieldList<DIM, std::vector<Vector> >&", "faceForce", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SVPHFacetedHydroBase)

