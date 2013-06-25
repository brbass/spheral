#-------------------------------------------------------------------------------
# ShadowConstantStrength
#
# Provides convenient constructors for the ConstantStrength model using the canned
# values in MaterialPropertiesLib.py.
#-------------------------------------------------------------------------------
from SpheralModules.Spheral.SolidMaterial import \
    ConstantStrength1d as RealConstantStrength1d, \
    ConstantStrength2d as RealConstantStrength2d, \
    ConstantStrength3d as RealConstantStrength3d
from SpheralModules.Spheral.Material import PhysicalConstants
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

#-------------------------------------------------------------------------------
# Define a string providing the help for building a ConstantStrength.
#-------------------------------------------------------------------------------
expectedUsageString = """
ConstantStrength can be constructed one of two ways:

1.  Using canned material values stored in our internal data base.  Expected arguments:
        materialName        : Label for the material in data base

2.  You can directly set the strength parameters explicitly, as
        mu0                 : shear modulus
        Y0                  : yield strength for plastic yielding
"""

#-------------------------------------------------------------------------------
# The base units for parameters in this file.
#-------------------------------------------------------------------------------
CGS = PhysicalConstants(0.01,    # Length in m
                        0.001,   # Mass in kg
                        1.0)     # Time in sec

#-------------------------------------------------------------------------------
# The generic factory function, where you pass in the dimension specific 
# ConstantStrength constructor.
# This one is for internal use only -- people will actually call the dimension
# specific front-ends at the end of this script.
#-------------------------------------------------------------------------------
def _ConstantStrengthFactory(*args, 
                              **kwargs):

    # The calling routine must provide the appropriate C++ constructor.
    CSConstructor = kwargs["CSConstructor"]

    # The arguments that need to be passed to this method.
    expectedArgs = ["materialName", "units"]
    optionalKwArgs = {"mu0" : None,
                      "Y0" : None}

    # What sort of information did the user pass in?
    if ("materialName" in kwargs or 
        len(args) > 0 and type(args[0]) is str):

        # It looks like the user is trying to use one of the libarary canned values.
        # Evaluate the arguments to the method.
        if len(args) > 0:
            if len(args) != len(expectedArgs):
                raise ValueError, expectedUsageString
            for i in xrange(len(expectedArgs)):
                exec("%s = args[i]" % expectedArgs[i])
            for arg in optionalKwArgs:
                exec("%s = optionalKwArgs['%s']" % (arg, arg))
        else:
            for arg in kwargs:
                if arg not in (expectedArgs + optionalKwArgs.keys() + ["CSConstructor"]):
                    raise ValueError, expectedUsageString
                exec("%s = kwargs['%s']" % (arg, arg))
            for arg in optionalKwArgs:
                if arg not in kwargs:
                    exec("%s = optionalKwArgs['%s']" % (arg, arg))

        # Check that the caller specified a valid material label.
        mat = materialName.lower()
        if mat not in SpheralMaterialPropertiesLib:
            raise ValueError, "You must specify one of %s" % str(SpheralMaterialPropertiesLib.keys())
        if ("mu0" not in SpheralMaterialPropertiesLib[mat] or
            "Y0"  not in SpheralMaterialPropertiesLib[mat]):
            raise ValueError, "The material %s does not provide strength paramters." % materialName

        # Extract the parameters for this material.
        if mu0 is None:
            mu0 = SpheralMaterialPropertiesLib[mat]["mu0"]
        if Y0 is None:
            Y0 = SpheralMaterialPropertiesLib[mat]["mu0"]
    
        # Figure out the conversions to the requested units.
        lconv = CGS.unitLengthMeters / units.unitLengthMeters
        mconv = CGS.unitMassKg / units.unitMassKg
        tconv = CGS.unitTimeSec / units.unitTimeSec
        rhoConv = mconv/(lconv*lconv*lconv)
        Pconv = mconv/(lconv*tconv*tconv)
        specificEconv = (lconv/tconv)**2

        # Build the arguments for constructing the ConstantStrength.
        passargs = [mu0 * Pconv,
                    Y0 * Pconv]
        passkwargs = {}

    else:

        # Just pass through the arguments.
        passargs = args
        passkwargs = kwargs
        del passkwargs["CSConstructor"]
    
    # Return the EOS.
    return CSConstructor(*tuple(passargs), **passkwargs)

#-------------------------------------------------------------------------------
# Create the dimension specific ConstantStrength factories.  These are the ones
# you actually use.
#-------------------------------------------------------------------------------
# 1D
def ConstantStrength1d(*args, **kwargs):
    kwargs["CSConstructor"] = RealConstantStrength1d
    return _ConstantStrengthFactory(*args, **kwargs)

# 2D
def ConstantStrength2d(*args, **kwargs):
    kwargs["CSConstructor"] = RealConstantStrength2d
    return _ConstantStrengthFactory(*args, **kwargs)

# 3D
def ConstantStrength3d(*args, **kwargs):
    kwargs["CSConstructor"] = RealConstantStrength3d
    return _ConstantStrengthFactory(*args, **kwargs)

#-------------------------------------------------------------------------------
# Override the doc strings to provide useful help.
#-------------------------------------------------------------------------------
ConstantStrength1d.__doc__ = expectedUsageString
ConstantStrength2d.__doc__ = expectedUsageString
ConstantStrength3d.__doc__ = expectedUsageString
