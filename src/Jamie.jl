module Jamie

using StaticArrays
using LinearAlgebra
using ForwardDiff
using DifferentialEquations

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                         EXPORT STATEMENTS                                        #
# ************************************************************************************************ #
# ************************************************************************************************ #

# Common
export name

# Util
export points_on_sphere

# Position-Velocity
export PositionVelocity
export position_magnitude, position_direction # position extended from Base.position
export velocity, velocity_magnitude, velocity_direction

# Dimensional Set
export DimensionalSet
export characteristic_mass, characteristic_length, characteristic_time
export characteristic_velocity, characteristic_acceleration
export dimensional_set

# Bodies
export mean_radius, mean_density, mean_eccentricity, mean_semimajor_axis
export gravitational_parameter, spice_identifier, equatorial_radius, polar_radius, ellipticity
export number_natural_satellites, solar_irradiance, total_mass, total_volume
export CelestialBody, parent_body, celestial_body, synodic_period
export NullShapeModel
export SphericalShapeModel
export NullPotential
export PointMassPotential

# Integration/Flow Map
export FlowMap
export dynamical_system, default_solver_options, default_parameters, solver_algorithm

# CRTBP
export CrtbpPrimary, CrtbpP1, CrtbpP2
export CrtbpSystem
export mass_ratio
export dimensionalize, nondimensionalize
export distance_from_primary
export pseudopotential, pseudopotential_gradient, pseudopotential_hessian
export jacobi_constant, jacobi_constant_gradient
export crtbp_eom, crtbp_jacobian
export equilibrium_solutions


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                          DIMENSIONAL SETS                                        #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    DimensionalSet(mass::M, length::L, time::T)

Structure defining a set of dimensional quantities for mass, length, and time.

See also:
[`characteristic_mass`](@ref), 
[`characteristic_time`](@ref),
[`characteristic_length`](@ref), 
[`characteristic_velocity`](@ref),
[`characteristic_acceleration`](@ref)
"""
struct DimensionalSet{M, L, T}
    mass::M
    length::L
    time::T
end

function DimensionalSet(;mass::M, length::L, time::T) where {M, L, T}
    DimensionalSet(mass, length, time)
end

"""
    characteristic_mass(::DimensionalSet)

Retrieve the characteristic mass of the dimensional set specified. 

See also: 
[`DimensionalSet`](@ref),
[`characteristic_time`](@ref),
[`characteristic_length`](@ref), 
[`characteristic_velocity`](@ref),
[`characteristic_acceleration`](@ref)
"""
characteristic_mass(ds::DimensionalSet) = ds.mass

"""
    characteristic_length(::DimensionalSet)

Retrieve the characteristic length of the dimensional set specified. 

See also: 
[`DimensionalSet`](@ref),
[`characteristic_mass`](@ref), 
[`characteristic_time`](@ref),
[`characteristic_velocity`](@ref),
[`characteristic_acceleration`](@ref)
"""
characteristic_length(ds::DimensionalSet) = ds.length

"""
    characteristic_time(::DimensionalSet)

Retrieve the characteristic time of the dimensional set specified. 

See also: 
[`DimensionalSet`](@ref),
[`characteristic_mass`](@ref), 
[`characteristic_length`](@ref),
[`characteristic_velocity`](@ref),
[`characteristic_acceleration`](@ref)
"""
characteristic_time(ds::DimensionalSet) = ds.time

"""
    characteristic_velocity(::DimensionalSet)

Retrieve the characteristic velocity of the dimensional set specified. 
Note, this is a derived quantity and equal to the characteristic length over the 
characteristic time.

See also: 
[`DimensionalSet`](@ref),
[`characteristic_mass`](@ref), 
[`characteristic_length`](@ref),
[`characteristic_time`](@ref),
[`characteristic_acceleration`](@ref)
"""
function characteristic_velocity(ds::DimensionalSet)
    len = characteristic_length(ds)
    t = characteristic_time(ds)
    len / t
end

"""
    characteristic_acceleration(::DimensionalSet)

Retrieve the characteristic acceleration of the dimensional set specified. 
Note, this is a derived quantity and equal to the characteristic length over the 
characteristic time squared.

See also: 
[`DimensionalSet`](@ref),
[`characteristic_mass`](@ref), 
[`characteristic_length`](@ref),
[`characteristic_time`](@ref),
[`characteristic_velocity`](@ref)
"""
function characteristic_acceleration(ds::DimensionalSet) 
    characteristic_velocity(ds) / characteristic_time(ds)
end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                           FILE INCLUDES                                          #
# ************************************************************************************************ #
# ************************************************************************************************ #
include("posvel.jl")
include("util.jl")
include("body/body.jl")
include("kepler.jl")
include("crtbp/crtbp.jl")
include("constants/constants.jl")
include("integration.jl")



end # module
