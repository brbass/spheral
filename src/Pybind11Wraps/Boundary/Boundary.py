#-------------------------------------------------------------------------------
# Boundary base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from BoundaryAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralBoundary")
class Boundary:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef typename Boundary<%(Dimension)s>::BoundaryNodes BoundaryNodes;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def setAllGhostNodes(self,
                         dataBase = "DataBase<%(Dimension)s>&"):
        "Recreate ghost nodes for this Boundary for all NodeLists in the DataBase."
        return "void"

    @PYB11virtual
    def setAllViolationNodes(self,
                             dataBase = "DataBase<%(Dimension)s>&"):
        "Select any nodes based in the FluidNodeLists in the given DataBase that are in violation of boundary condition."
        return "void"

    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<%(Dimension)s, int>&",
                       old2newIndexMap = "FieldList<%(Dimension)s, int>&",
                       numNodesRemoved = "std::vector<int>&"):
        "Use a set of flags to cull out inactive ghost nodes."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary20(self,
                             field = "Field<%(Dimension)s, std::vector<Scalar>>&"):
        "We provide default copies for arrays of values, but descendants can override these."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary21(self,
                             field = "Field<%(Dimension)s, std::vector<Vector>>&"):
        "We provide default copies for arrays of values, but descendants can override these."
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self):
        "Some boundaries need to know when a problem is starting up and all the physics packages have been initialized."
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeGhostBoundary(self):
        "Provide an optional hook that is to be called when all ghost boundaries are to have been set."
        return "void"

    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<%(Dimension)s>&"):
        "Overridable hook for clearing out the boundary condition."
        return "void"

    @PYB11virtual
    @PYB11const
    def clip(self,
             xmin = "Vector&",
             xmax = "Vector&"):
        "Optionally the boundary can clip an input box range. Defaults to no-op."
        return "void"

    @PYB11virtual
    @PYB11const
    def allowGhostCulling(self):
        "Optionally opt-out of ghost node culling."
        return "bool"

    @PYB11virtual
    @PYB11const
    def meshGhostNodes(self):
        "Some boundaries have ghosts we should exclude from tessellations. Provide a hook to note such cases."
        return "bool"

    @PYB11virtual
    def addNodeList(self, nodeList="NodeList<%(Dimension)s>&"):
        return "void"
                    
    #...........................................................................
    # Methods
    @PYB11const
    def haveNodeList(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Check if we have entries for the given NodeList."
        return "bool"

    @PYB11returnpolicy("reference_internal")
    def accessBoundaryNodes(self, nodeList="NodeList<%(Dimension)s>&"):
        return "BoundaryNodes&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def controlNodes(self, nodeList="NodeList<%(Dimension)s>&"):
        "Control nodes for a given NodeList"
        return "const std::vector<int>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def ghostNodes(self, nodeList="NodeList<%(Dimension)s>&"):
        "Ghost nodes for a given NodeList"
        return "const std::vector<int>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def violationNodes(self, nodeList="NodeList<%(Dimension)s>&"):
        "Violation nodes for a given NodeList"
        return "const std::vector<int>&"

    #...........................................................................
    # applyFieldListGhostBoundary/enforceFieldListBoundary
    @PYB11template("Value")
    @PYB11const
    def applyFieldListGhostBoundary(self, fieldList="FieldList<%(Dimension)s, %(Value)s>&"):
        "Apply ghost boundary to all Fields in FieldList"
        return "void"

    @PYB11template("Value")
    @PYB11const
    def enforceFieldListBoundary(self, fieldList="FieldList<%(Dimension)s, %(Value)s>&"):
        "Enforce boundary on all Fields in FieldList"
        return "void"

    for T in ("int", "Scalar", "Vector", "Tensor", "SymTensor", "ThirdRankTensor", "FourthRankTensor", "FifthRankTensor"):
        exec('''
aflgb%(T)s = PYB11TemplateMethod(applyFieldListGhostBoundary, template_parameters="%(T)s", pyname="applyFieldListGhostBoundary")
eflgb%(T)s = PYB11TemplateMethod(enforceFieldListBoundary, template_parameters="%(T)s", pyname="enforceFieldListBoundary")
''' % {"T" : T})

    #...........................................................................
    # Attributes
    numGhostNodes = PYB11property("int")

    #---------------------------------------------------------------------------
    # BoundaryNodes
    #---------------------------------------------------------------------------
    class BoundaryNodes:

        controlNodes = PYB11readwrite()
        ghostNodes = PYB11readwrite()
        violationNodes = PYB11readwrite()

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, Boundary, pure_virtual=True)
