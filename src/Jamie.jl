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
export pv_position, pv_position_mag, pv_position_unit
export pv_velocity, pv_velocity_mag, pv_velocity_unit

# Dimensional Set
export DimensionalSet
export characteristic_mass, characteristic_length, characteristic_time
export characteristic_velocity, characteristic_acceleration

# Bodies
export name_string
export gravitational_parameter, gravitational_potential, gravitational_acceleration
export mean_radius
export CelestialBody, string, spice_identifier, potential_model, shape_model
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
#                                           INTEGRATION PROBLEMS                                   #
# ************************************************************************************************ #
# ************************************************************************************************ #
const DEFAULT_SOLVE_OPTS = Dict(
    :abstol => 1e-12,
    :reltol => 1e-12,
)

struct Propagation{SYS, A, S, P}
    system::SYS
    algorithm::A
    solve_opts::S
    p::P
end

function Propagation(sys; algorithm=Vern9(), solve_opts=DEFAULT_SOLVE_OPTS, params=nothing)
    Propagation(
        sys,
        algorithm,
        solve_opts,
        params
    )
end

equations_of_motion(prop::Propagation) = prop.system
default_parameters(prop::Propagation) = prop.p
solver_options(prop::Propagation) = prop.solve_opts
algorithm(prop::Propagation) = prop.algorithm


function (prop::Propagation)(q0, tspan; p::P=nothing, solve_opts...) where {P}
    if P == Nothing
        params = default_parameters(prop)
    else
        params = p
    end

    prob = ODEProblem(
        equations_of_motion(prop),
        q0,
        tspan,
        params
    )
    solve(
        prob,
        algorithm(prop);
        solver_options(prop)...,
        solve_opts...
    )
end

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                           FILE INCLUDES                                          #
# ************************************************************************************************ #
# ************************************************************************************************ #
include("util.jl")
include("body/body.jl")
include("kepler.jl")
include("crtbp.jl")
include("constants.jl")



end # module
