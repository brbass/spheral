#-------------------------------------------------------------------------------
# GridCellIndex template
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class GridCellIndex:

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self, xIndex="int"):
        "Construct with (x)"

    def pyinit2(self, xIndex="int", yIndex="int"):
        "Construct with (x,y)"

    def pyinit3(self, xIndex="int", yIndex="int", zIndex="int"):
        "Construct with (x,y,z)"

    def pyinit4(self, rhs="const GridCellIndex<%(Dimension)s>&"):
        "Copy constructor"

    #...........................................................................
    # Methods
    @PYB11pycppname("setIndices")
    def setIndices1(self, xIndex="int"):
        "Set the x index"
        return "void"

    @PYB11pycppname("setIndices")
    def setIndices2(self, xIndex="int", yIndex="int"):
        "Set the (x,y) indices"
        return "void"

    @PYB11pycppname("setIndices")
    def setIndices3(self, xIndex="int", yIndex="int", zIndex="int"):
        "Set the (x,y,z) indices"
        return "void"

    def dot(self):
        "Dot product with another GridCellIndex"
        return

    #...........................................................................
    # Sequence methods
    @PYB11implementation("[](const GridCellIndex<%(Dimension)s>& self) { return %(Dimension)s::nDim; }")
    def __len__(self):
        "The size (in number of coordinates) of the GridCellIndex."

    @PYB11implementation("[](const GridCellIndex<%(Dimension)s> &s, int i) { if (i >= %(Dimension)s::nDim) throw py::index_error(); return s(i); }") 
    def __getitem__(self):
        "Python indexing to get a coordinate."

    @PYB11implementation("[](GridCellIndex<%(Dimension)s> &s, int i, float v) { if (i >= %(Dimension)s::nDim) throw py::index_error(); s(i) = v; }") 
    def __setitem__(self):
        "Python indexing to set a coordinate."

    @PYB11implementation("[](const GridCellIndex<%(Dimension)s> &s) { return py::make_iterator(std::begin(s), std::end(s)); }")
    def __iter__(self):
        "Python iteration through a GridCellIndex."

    @PYB11const
    def __call__(self, i="int"):
        "Index for a coordinate using parens."
        return "int"
