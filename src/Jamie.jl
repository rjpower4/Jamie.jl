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
export name, mean_radius, mean_density, mean_eccentricity, mean_semimajor_axis
export gravitational_parameter, spice_identifier, equatorial_radius, polar_radius, ellipticity
export number_natural_satellites, solar_irradiance, total_mass, total_volume
export CelestialBody, parent_body, celestial_body, synodic_period
export NullShapeModel
export SphericalShapeModel
export NullPotential
export PointMassPotential

# CRTBP
export CrtbpPrimary, CrtbpP1, CrtbpP2
export CrtbpSystem, mass_ratio 
export dimensionalize, nondimensionalize
export distance_from_primary
export pseudopotential, pseudopotential_gradient, pseudopotential_hessian
export jacobi_constant, jacobi_constant_gradient
export crtbp_eom, crtbp_jacobian
export equilibrium_solutions


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                         POSITION-VELOCITY                                        #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    PositionVelocity{T} <: FieldVector{6, T}

Struct representing the position-velocity vector.
`T` designates the types of component in the vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
struct PositionVelocity{T} <: FieldVector{6, T}
    x::T
    y::T
    z::T
    vx::T
    vy::T
    vz::T
end

"""
    position(::PositionVelocity)

Get the position components of the position-velocity vector.

See also: [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
Base.position(pv::PositionVelocity) = @SVector [ pv.x, pv.y, pv.z ]

"""
    position_magnitude(::PositionVelocity)

Return the magnitude of the position components of the position-velocity vector.
This is equivalent to the distance from the origin of the position vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
position_magnitude(pv::PositionVelocity) = (norm ∘ position)(pv)

"""
    position_direction(::PositionVelocity)

Return the unit vector parallel to the position sub-vector in the specified position-velocity vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
position_direction(pv::PositionVelocity) = position(pv) ./ position_magnitude(pv)

"""
    velocity(::PositionVelocity)

Get the velocity components of the position-velocity vector.

See also: [`position`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
velocity(pv::PositionVelocity) = @SVector [ pv.vx, pv.vy, pv.vz ]

"""
    velocity_magnitude(::PositionVelocity)

Return the magnitude of the velocity components of the position-velocity vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref),
[`velocity_direction`](@ref)
"""
velocity_magnitude(pv::PositionVelocity) = (norm ∘ velocity)(pv)

"""
    velocity_direction(::PositionVelocity)

Return the unit vector parallel to the velocity sub-vector in the specified position-velocity vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
"""
velocity_direction(pv::PositionVelocity) = velocity(pv) ./ velocity_magnitude(pv)


"""
    dualize(pv::PositionVelocity)

Build a dual position-velocity to enable evaluation of partial derivatives.
"""
function ForwardDiff.dualize(pv::PositionVelocity{T}) where {T}
    PositionVelocity(
        ForwardDiff.dualize(PositionVelocity{T}, pv)
    )
end

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
include("util.jl")
include("body/body.jl")
include("kepler.jl")
include("crtbp/crtbp.jl")
include("constants/constants.jl")
include("integration.jl")



end # module
