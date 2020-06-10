module Jamie

using StaticArrays
using LinearAlgebra

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                         EXPORT STATEMENTS                                        #
# ************************************************************************************************ #
# ************************************************************************************************ #

# Util
export points_on_sphere

# Position-Velocity
export PositionVelocity
export pv_position, pv_position_mag, pv_position_unit
export pv_velocity, pv_velocity_mag, pv_velocity_unit

# Bodies
export gravitational_parameter, mean_radius
export CelestialBody, name_string, spice_identifier, potential_model, shape_model
export NullShapeModel
export SphericalShapeModel
export NullPotential
export PointMassPotential

# CRTBP
export CrtbpPrimary, CrtbpP1, CrtbpP2
export CrtbpSystem, mass_ratio, characteristic_mass, characteristic_length, characteristic_time
export characteristic_velocity
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

See also: [`pv_position`](@ref), [`pv_velocity`](@ref),
[`pv_position_mag`](@ref), [`pv_position_unit`](@ref), [`pv_velocity_mag`](@ref),
[`pv_velocity_unit`](@ref)
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
    pv_position(::PositionVelocity)

Get the position components of the position-velocity vector.

See also: [`pv_velocity`](@ref),
[`pv_position_mag`](@ref), [`pv_position_unit`](@ref), [`pv_velocity_mag`](@ref),
[`pv_velocity_unit`](@ref)
"""
pv_position(pv::PositionVelocity) = @SVector [ pv.x, pv.y, pv.z ]

"""
    pv_position_mag(::PositionVelocity)

Return the magnitude of the position components of the position-velocity vector.
This is equivalent to the distance from the origin of the position vector.

See also: [`pv_position`](@ref), [`pv_velocity`](@ref),
[`pv_position_unit`](@ref), [`pv_velocity_mag`](@ref),
[`pv_velocity_unit`](@ref)
"""
pv_position_mag(pv::PositionVelocity) = (norm ∘ pv_position)(pv)

"""
    pv_position_unit(::PositionVelocity)

Return the unit vector parallel to the position sub-vector in the specified position-velocity vector.

See also: [`pv_position`](@ref), [`pv_velocity`](@ref),
[`pv_position_mag`](@ref), [`pv_velocity_mag`](@ref),
[`pv_velocity_unit`](@ref)
"""
pv_position_unit(pv::PositionVelocity) = pv_position(pv) ./ pv_position_mag(pv)

"""
    pv_velocity(::PositionVelocity)

Get the velocity components of the position-velocity vector.

See also: [`pv_position`](@ref),
[`pv_position_mag`](@ref), [`pv_position_unit`](@ref), [`pv_velocity_mag`](@ref),
[`pv_velocity_unit`](@ref)
"""
pv_velocity(pv::PositionVelocity) = @SVector [ pv.vx, pv.vy, pv.vz ]

"""
    pv_velocity_mag(::PositionVelocity)

Return the magnitude of the velocity components of the position-velocity vector.

See also: [`pv_position`](@ref), [`pv_velocity`](@ref),
[`pv_position_mag`](@ref), [`pv_position_unit`](@ref),
[`pv_velocity_unit`](@ref)
"""
pv_velocity_mag(pv::PositionVelocity) = (norm ∘ pv_velocity)(pv)

"""
    pv_velocity_unit(::PositionVelocity)

Return the unit vector parallel to the velocity sub-vector in the specified position-velocity vector.

See also: [`pv_position`](@ref), [`pv_velocity`](@ref),
[`pv_position_mag`](@ref), [`pv_position_unit`](@ref), [`pv_velocity_mag`](@ref),
"""
pv_velocity_unit(pv::PositionVelocity) = pv_velocity(pv) ./ pv_velocity_mag(pv)


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                           FILE INCLUDES                                          #
# ************************************************************************************************ #
# ************************************************************************************************ #
include("util.jl")
include("body.jl")
include("crtbp.jl")



end # module
