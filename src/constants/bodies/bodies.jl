struct Sun <: CelestialBody end
name(::Sun) = "SUN"
spice_identifier(::Sun) = 399
mean_radius(::Sun) = 695700.0 # km
ellipticity(::Sun) = 0.00005 
mean_density(::Sun) = 1408.0 # kg / m^3
gravitational_parameter(::Sun) = 1.32712e11 # km^3 / s^2
total_mass(::Sun) = 1.9885e30 # kg
total_volume(::Sun) = 1.412e18 # km^3


include("earth-system.jl")