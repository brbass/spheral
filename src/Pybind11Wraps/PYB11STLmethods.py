#-------------------------------------------------------------------------------
# PYB11STLmethods
#
# Thin wrappers to generate the pybind11 STL functions.
#-------------------------------------------------------------------------------
import inspect
from PYB11utils import *

#-------------------------------------------------------------------------------
# std::vector
#-------------------------------------------------------------------------------
class PYB11_bind_vector:

    def __init__(self,
                 element,            # template element type for vector
                 opaque = False):    # should we make the type opaque
        self.element = element
        self.opaque = opaque
        return

    def preamble(self, modobj, ss, name):
        if self.opaque:
            if "," in self.element:
                ss('''{
  typedef %(element)s %(safe_element)s; 
  PYBIND11_MAKE_OPAQUE(std::vector<%(safe_element)s>);
}\n''' % {"element": self.element,
          "safe_element" : PYB11mangle(self.element)})
            else:
                ss('PYBIND11_MAKE_OPAQUE(std::vector<' + self.element + '>);\n')
        return

    def __call__(self, modobj, ss, name):
        ss('py::bind_vector<std::vector<' + self.element + '>>(m, "' + name + '");\n')
        return

#-------------------------------------------------------------------------------
# std::map
#-------------------------------------------------------------------------------
class PYB11_bind_map:

    def __init__(self,
                 key,                # template key type
                 value,              # template value type
                 opaque = False):    # should we make the container opaque
        self.key = key
        self.value = value
        self.opaque = opaque
        return

    def preamble(self, modobj, ss, name):
        if self.opaque:
            cppname = "std::map<" + self.key + ", " + self.value + ">"
            manglename = PYB11mangle(cppname)
            ss('''{
  typedef %(cppname)s %(manglename)s;
  PYBIND11_MAKE_OPAQUE(%(manglename)s);
}
''' % {"cppname" : cppname,
       "manglename" : manglename})
        return

    def __call__(self, modobj, ss, name):
        ss('py::bind_map<std::map<' + self.key + ', ' + self.value + '>>(m, "' + name + '");\n')
        return

#-------------------------------------------------------------------------------
# PYB11STLobjs
#
# Get the STL objects to bind from a module
#-------------------------------------------------------------------------------
def PYB11STLobjs(modobj):
    return [(name, obj) for (name, obj) in inspect.getmembers(modobj)
            if name[:5] != "PYB11" and
            (isinstance(obj, PYB11_bind_vector) or
             isinstance(obj, PYB11_bind_map))]

#-------------------------------------------------------------------------------
# PYB11generateModuleSTL
#
# Bind the STL containers in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleSTL(modobj, ss):
    stuff = PYB11STLobjs(modobj)
    for (name, obj) in stuff:
        ss("  ")
        obj(modobj, ss, name)
    ss("\n")
    return

