"""
    GRAVITATIONAL_CONSTANT

The standard Newtonian gravitational constant in km^3 / (kg s^2).
The value is 6.674e-20.
"""
const GRAVITATIONAL_CONSTANT = 6.6743e-20 # km^3 / (kg s^2)

"""
    SPEED_OF_LIGHT

The speed of light in a vacuum in km/s as defined by NIST.
"""
const SPEED_OF_LIGHT = 299792.458 # km / s

include("integration.jl")
include("bodies/bodies.jl")