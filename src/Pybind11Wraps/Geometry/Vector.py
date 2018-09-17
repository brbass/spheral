from PYB11Decorators import *
from PYB11property import *
from PYB11class import *

#-------------------------------------------------------------------------------
# Vector template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class Vector:
    "Spheral Geometric Vector class"

    # Static attributes
    @PYB11static
    @PYB11readonly
    def nDimensions(self):
        "Number of dimensions"

    @PYB11static
    @PYB11readonly
    def numElements(self):
        "Number of elements stored in the type."

    @PYB11static
    @PYB11readonly
    def zero(self):
        "The zero value equivalent."

    @PYB11static
    @PYB11readonly
    def one(self):
        "The unit value equivalent."

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Dim<%(ndim)s>::Vector"):
        "Copy constructor"

    def pyinit2(self,
                x = "double",
                y = ("double", "0.0"),
                z = ("double", "0.0")):
        "Construct with element values."

    # Methods
    @PYB11const
    def dot(self):
        "Dot (inner) product with a Vector."

    @PYB11const
    def cross(self):
        "Cross product with a Vector."

    @PYB11const
    def dyad(self):
        "Dyadic (outer) product with a Vector."

    @PYB11const
    def selfdyad(self):
        "Dyadic (outer) product with ourself."

    @PYB11const
    def unitVector(self):
        "Unit vector in the direction of this one."

    @PYB11const
    def magnitude(self):
        "The magnitude of the Vector."

    @PYB11const
    def magnitude2(self):
        "The square of the magnitude of the Vector."
        return

    @PYB11const
    def minElement(self):
        "Minimum (x,y,z) in the Vector."

    @PYB11const
    def maxElement(self):
        "Maximum (x,y,z) in the Vector."

    @PYB11const
    def maxAbsElement(self):
        "Maximum absolute element (|x|,|y|,|z|) in the Vector."

    @PYB11const
    def sumElements(self):
        "Sum of the elements (x+y+z) in the Vector."

    # Operators
    def __neg__(self):
        return
    def __add__(self):
        return
    def __sub__(self):
        return
    def __mul__(self):
        return
    def __iadd__(self):
        return
    def __isub__(self):
        return
    def __mul__(self, rhs="float()"):
        return
    def __rmul__(self, rhs="float()"):
        return
    def __div__(self, rhs="float()"):
        return
    def __imul__(self, rhs="float()"):
        return
    def __idiv__(self, rhs="float()"):
        return

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::Vector& self) { return Dim<%(ndim)s>::nDim; }")
    def __len__(self):
        "The size (in number of coordinates) of the Vector."

    @PYB11implementation("[](const Dim<%(ndim)s>::Vector &s, size_t i) { if (i >= Dim<%(ndim)s>::Vector::numElements) throw py::index_error(); return s[i]; }") 
    def __getitem__(self):
        "Python indexing to get a coordinate."

    @PYB11implementation("[](Dim<%(ndim)s>::Vector &s, size_t i, float v) { if (i >= Dim<%(ndim)s>::Vector::numElements) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set a coordinate."

    @PYB11implementation("[](const Dim<%(ndim)s>::Vector &s) { return py::make_iterator(s.begin(), s.end()); }")
    def __iter__(self):
        "Python iteration through a Vector."

    @PYB11const
    def __call__(self, i="Dim<%(ndim)s>::Vector::size_type"):
        "Index for a coordinate using parens."
        return "double"

    # Comparison
    @PYB11const
    def compare(self,
                rhs = "const Dim<%(ndim)s>::Vector&"):
        "Compare (-1,0,1) with a Vector."
        return "int"

    @PYB11const
    @PYB11pyname("compare")
    @PYB11cppname("compare")
    def compare1(self,
                 rhs = "const double"):
        "Compare (-1,0,1) with a double."
        return "int"

    def __eq__(self):
        return
    def __ne__(self):
        return
    def __lt__(self):
        return
    def __gt__(self):
        return
    def __le__(self):
        return
    def __ge__(self):
        return

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::Vector& self) {
  std::string result = "Vector" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

    # x
    @PYB11cppname("x")
    @PYB11const
    @PYB11ignore
    def getx(self):
        return "double"

    @PYB11cppname("x")
    @PYB11ignore
    def setx(self, val="double"):
        return "void"

    # y
    @PYB11cppname("y")
    @PYB11const
    @PYB11ignore
    def gety(self):
        return "double"

    @PYB11cppname("y")
    @PYB11ignore
    def sety(self, val="double"):
        return "void"

    # z
    @PYB11cppname("z")
    @PYB11const
    @PYB11ignore
    def getz(self):
        return "double"

    @PYB11cppname("z")
    @PYB11ignore
    def setz(self, val="double"):
        return "void"

    # Properties
    x = property(getx, setx, doc="The x coordinate.")
    y = property(gety, sety, doc="The y coordinate.")
    z = property(getz, setz, doc="The z coordinate.")

#-------------------------------------------------------------------------------
# Vector instantiations.
#-------------------------------------------------------------------------------
Vector1d = PYB11TemplateClass(Vector,
                              template_parameters = ("1"),
                              cppname = "Dim<1>::Vector",
                              pyname = "Vector1d",
                              docext = " (1D).")
Vector2d = PYB11TemplateClass(Vector,
                              template_parameters = ("2"),
                              cppname = "Dim<2>::Vector",
                              pyname = "Vector2d",
                              docext = " (2D).")
Vector3d = PYB11TemplateClass(Vector,
                              template_parameters = ("3"),
                              cppname = "Dim<3>::Vector",
                              pyname = "Vector3d",
                              docext = " (3D).")

