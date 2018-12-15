#-------------------------------------------------------------------------------
# FileIO abstract class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIOAbstractMethods import *
from FileIOTemplateMethods import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11module("SpheralFileIO")
class FileIO:
    "Abstract base class for FileIO objects"

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Construct with a given file name and access"

    #...........................................................................
    # Abstract methods
    @PYB11pure_virtual
    def open(self,
             fileName = "const std::string",
             access = "AccessType"):
        "Open a file for IO"
        return "void"

    @PYB11pure_virtual
    def close(self):
        "Close the current file we're pointing at"
        return "void"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def write_unsigned_int(self,
                           value = "const unsigned int",
                           pathName = "const std::string"):
        "Write an unsigned int"
        return "void"

    @PYB11virtual
    def write_int(self,
                  value = "const int",
                  pathName = "const std::string"):
        "Write an int"
        return "void"

    @PYB11virtual
    def write_bool(self,
                   value = "const bool",
                   pathName = "const std::string"):
        "Write a bool"
        return "void"

    @PYB11virtual
    def write_double(self,
                     value = "const double",
                     pathName = "const std::string"):
        "Write a double"
        return "void"

    @PYB11virtual
    def write_string(self,
                     value = "const std::string",
                     pathName = "const std::string"):
        "Write a std::string"
        return "void"

    @PYB11virtual
    @PYB11const
    def read_unsigned_int(self,
                          pathName = "const std::string"):
        "Read an unsigned int"
        return "unsigned int"

    @PYB11virtual
    @PYB11const
    def read_int(self,
                 pathName = "const std::string"):
        "Read an int"
        return "int"

    @PYB11virtual
    @PYB11const
    def read_bool(self,
                  pathName = "const std::string"):
        "Read a bool"
        return "bool"

    @PYB11virtual
    @PYB11const
    def read_double(self,
                    pathName = "const std::string"):
        "Read a double"
        return "double"

    @PYB11virtual
    @PYB11const
    def read_string(self,
                    pathName = "const std::string"):
        "Read a std::string"
        return "std::string"

    #...........................................................................
    # Methods
    for ndim in xrange(1,4):   # all three always required
        exec('''
@PYB11pycppname("write")
def writePlane%(ndim)i(self,
                       value = "const GeomPlane<Dim<%(ndim)i>>&",
                       pathName = "const std::string"):
    "Write a Plane%(ndim)id"
    return "void"

@PYB11pycppname("read")
@PYB11const
def readPlane%(ndim)i(self,
                      value = "GeomPlane<Dim<%(ndim)i>>&",
                      pathName = "const std::string"):
    "Read a Plane%(ndim)id"
    return "void"
''' % {"ndim" : ndim})

    @PYB11const
    def splitPathComponents(self, pathName="const std::string&"):
        "A helper function to split a string up into substrings delimited by '/'."
        return "std::vector<std::string>"

    @PYB11const
    def groupName(self, pathName="const std::string"):
        "Return the group (directory) component of a path."
        return "std::string"

    @PYB11const
    def variableName(self, pathName="const std::string"):
        "Return the variable component of a path."
        return "std::string"

    @PYB11implementation("[](FileIO& self, py::handle thing, py::handle path) { self.writeObject(thing.ptr(), path.ptr()); }")
    def writeObject(self,
                    thing = "py::handle",
                    path = "py::handle"):
        "Handle a generic python object through serialization"
        return "void"

    @PYB11returnpolicy("take_ownership")
    @PYB11const
    @PYB11implementation("[](FileIO& self, py::handle path) { return py::handle(self.readObject(path.ptr())); }")
    def readObject(self,
                   path = "py::handle"):
        "Return a generic python object from deserialization."
        return "py::handle"

    #...........................................................................
    # Properties
    fileName = PYB11property("const std::string&", "fileName", doc="The current file name")
    access = PYB11property("AccessType", "access", doc="The access type of the currently open file")
    fileOpen = PYB11property("bool", "fileOpen", doc="Is the file currently open?")

#-------------------------------------------------------------------------------
# Inject the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, FileIO, virtual=False, pure_virtual=True)
PYB11inject(FileIOTemplateMethods, FileIO)