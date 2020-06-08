# ************************************************************************************************ #
# ************************************************************************************************ #
#                                       PRIMARY DEFINITIONS                                        #
# ************************************************************************************************ #
# ************************************************************************************************ #
abstract type CrtbpPrimary end

"""
    CrtbpP1

Empty struct identifying the first (usually larger) primary in the CRTBP.

See also: [`CrtbpP2`](@ref)
"""
struct CrtbpP1 <: CrtbpPrimary end

"""
    CrtbpP2

Empty struct identifying the second (usually smaller) primary in the CRTBP.

See also: [`CrtbpP1`](@ref)
"""
struct CrtbpP2 <: CrtbpPrimary end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                     CRTBP SYSTEM DEFINITIONS                                     #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    CrtbpSystem{T, D <: Union{T, Nothing}}(μ::T, name::String, char_mass::D, char_length::D, char_time::D)

Struct defining the parameters for the Crtbp System. Contains the mass ratio parameter, ``\\mu``,
as well as the system name and characteristic quantities.

See also: [`mass_ratio`](@ref), [`characteristic_mass`](@ref), [`characteristic_time`](@ref),
[`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
struct CrtbpSystem{T, D <: Union{T, Nothing}}
    μ::T
    name::String
    char_mass::D
    char_length::D
    char_time::D

    function CrtbpSystem(μ::T, name::String, char_mass::D, char_length::D, char_time::D) where {T,D}
        if μ <= 0.0
            throw(DomainError(μ, "Non-positive μ value"))
        end
        if D != Nothing
            if char_mass <= 0
                throw(DomainError(char_mass, "Non-Positive characteristic mass"))
            elseif char_length<= 0
                throw(DomainError(char_length, "Non-Positive characteristic length"))
            elseif char_time <= 0
                throw(DomainError(char_time, "Non-Positive characteristic time"))
            end
        end
        new{T, D}(μ, name, char_mass, char_length, char_time)
    end
end

name_string(s::CrtbpSystem) = s.name

"""
    CrtbpSystem(μ::T; [name, char_mass, char_length, char_time])

Construct a new `CrtbpSystem` with the μ value and the optional keywor arguments.

# Arguments
- `μ::T`: The mass ratio value of the system
- `name::String="UNNAMMED SYSTEM"`: name for the specified system
- `char_mass::Union{T, Nothing}=nothing`: characteristic mass for the system
- `char_length::Union{T, Nothing}=nothing`: characteristic length for the system
- `char_time::Union{T, Nothing}=nothing`: characteristic time for the system
"""
function CrtbpSystem(μ::T; name="UNNAMED SYSTEM",
                     char_mass::D=nothing,
                     char_length::D=nothing,
                     char_time::D=nothing) where {T, D}
    CrtbpSystem(μ, name, char_mass, char_length, char_time)
end


"""
    mass_ratio(::CrtbpSystem)

Retrieve the mass ratio, μ, of the specified system.

See also: [`CrtbpSystem`](@ref), [`characteristic_mass`](@ref), [`characteristic_time`](@ref),
[`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
mass_ratio(s::CrtbpSystem) = s.μ

"""
    characteristic_mass(::CrtbpSystem)

Retrieve the characteristic mass of the CRTBP system specified. This mass is usually defined as
the sum of the masses of the two primaries.

See also: [`CrtbpSystem`](@ref), [`mass_ratio`](@ref), [`characteristic_time`](@ref),
[`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
characteristic_mass(s::CrtbpSystem{T,T}) where {T} = s.char_mass

"""
    characteristic_length(::CrtbpSystem)

Retrieve the characteristic length of the CRTBP system specified. This length is usually defined as
the semi-major axis of the circular orbit of P2 in the CRTBP.

See also: [`CrtbpSystem`](@ref), [`mass_ratio`](@ref), [`characteristic_mass`](@ref),
[`characteristic_time`](@ref), [`characteristic_velocity`](@ref)
"""
characteristic_length(s::CrtbpSystem{T, T}) where {T} = s.char_length

"""
    characteristic_length(::CrtbpSystem)

Retrieve the characteristic time of the CRTBP system specified. This time is usually defined as the
time such that the universal gravitational constant in the non-dimensional system is unity.

See also: [`CrtbpSystem`](@ref), [`mass_ratio`](@ref), [`characteristic_mass`](@ref),
[`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
characteristic_time(s::CrtbpSystem{T, T}) where {T} = s.char_time

"""
    characteristic_velocity(::CrtbpSystem)

Retrieve the characteristic velocityof the CRTBP system specified. This is equivalent to the
characteristic length divided by the characteristic time.

See also: [`CrtbpSystem`](@ref), [`mass_ratio`](@ref), [`characteristic_mass`](@ref),
[`characteristic_length`](@ref)
"""
function characteristic_velocity(s::CrtbpSystem{T, T}) where {T}
    characteristic_length(s) / characteristic_time(s)
end

"""
    dimensionalize(::CrtbpSystem{T,T}, ::PositionVelocity)

Dimensionalize the position and velocity components using the characteristic quantities associated
with the CRTBP system.

See also: [`nondimensionalize`](@ref)
[`CrtbpSystem`](@ref), [`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
function dimensionalize(sys::CrtbpSystem{T1,T1}, pv::PositionVelocity{T2}) where {T1, T2}
    cl = characteristic_length(sys)
    cv = characteristic_velocity(sys)
    scale_vec = @SVector [cl, cl, cl, cv, cv, cv]
    pv .* scale_vec
end

"""
    nondimensionalize(::CrtbpSystem{T,T}, ::PositionVelocity)

Nondimensionalize the position and velocity components using the characteristic quantities associated
with the CRTBP system.

See also: [`dimensionalize`](@ref)
[`CrtbpSystem`](@ref), [`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
function nondimensionalize(sys::CrtbpSystem{T1,T1}, pv::PositionVelocity{T2}) where {T1, T2}
    cl = characteristic_length(sys)
    cv = characteristic_velocity(sys)
    scale_vec = @SVector [cl, cl, cl, cv, cv, cv]
    pv ./ scale_vec
end

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                   PRIMARY DISTANCE FUNCTIONS                                     #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    distance_from_primary(::CrtbpSystem, ::PositionVelocity)

Determine the scalar distance from the primaries in the CRTBP.
Values are returned as a named tuple with fields `P1` and `P2`.

See also: [`CrtbpSystem`](@ref), [`CrtbpP1`](@ref), [`CrtbpP2`](@ref)
"""
function distance_from_primary(sys::CrtbpSystem, pv::PositionVelocity)
    (
        P1=distance_from_primary(sys, pv, CrtbpP1()),
        P2=distance_from_primary(sys, pv, CrtbpP2())
    )
end

"""
    distance_from_primary(::CrtbpSystem, ::PositionVelocity, ::CrtbpP1)

Determine the scalar distance from the first primary in the CRTBP.
"""
function distance_from_primary(sys::CrtbpSystem, pv::PositionVelocity, ::CrtbpP1)
    xx = pv[1] + mass_ratio(sys)
    sqrt(xx^2 + pv[2]^2 + pv[3]^2)
end

"""
    distance_from_primary(::CrtbpSystem, ::PositionVelocity, ::CrtbpP2)

Determine the scalar distance from the second primary in the CRTBP.
"""
function distance_from_primary(sys::CrtbpSystem, pv::PositionVelocity, ::CrtbpP2)
    xx = pv[1] - 1.0 + mass_ratio(sys)
    sqrt(xx^2 + pv[2]^2 + pv[3]^2)
end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                         PSEUDOPOTENTIAL                                          #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    pseudopotential(::CrtbpSystem, ::PositionVelocity)

Calculate the CRTBP pseudopotential of the specified position-velocity vector.

See also: [`pseudopotential_gradient`](@ref), [`pseudopotential_hessian`](@ref)
"""
function pseudopotential(sys::CrtbpSystem, pv::PositionVelocity)
    μ = mass_ratio(sys)
    (r13, r23) = distance_from_primary(sys, pv)
    (1.0/2.0) * (pv[1]^2 + pv[2]^2) + (1.0 - μ) / r13 + μ / r23
end

"""
    pseudopotential_gradient(::CrtbpSystem, ::PositionVelocity)

Calculate the partials of the CRTBP pseudopotential with respect to all 6 state components in the
position-velocity vector.

See also: [`pseudopotential`](@ref), [`pseudopotential_hessian`](@ref)
"""
function pseudopotential_gradient(sys::CrtbpSystem, pv::PositionVelocity)
    μ = mass_ratio(sys)
    (r13, r23) = distance_from_primary(sys, pv)
    r13_3 = r13^3
    r23_3 = r23^3
    omm = 1.0 - μ
    @SVector [
        pv[1] - (omm) * (pv[1] + μ) / r13_3 - μ * (pv[1] - omm) / r23_3,
        pv[2] - (omm) * pv[2] / r13_3 - μ * pv[2] / r23_3,
        -(omm) * pv[3] / r13_3 - μ * pv[3] / r23_3,
        0.0,
        0.0,
        0.0
    ]
end

"""
    pseudopotential_hessian(::CrtbpSystem, ::PositionVelocity)

Calculate the second order partials of the CRTBP pseudopotential with respect to all 6 state
components in the position-velocity vector.

See also: [`pseudopotential`](@ref), [`pseudopotential_gradient`](@ref)
"""
function pseudopotential_hessian(sys::CrtbpSystem, pv::PositionVelocity)
    μ = mass_ratio(sys)

    (r13, r23) = distance_from_primary(sys, pv)
    r13_3 = r13^3
    r23_3 = r23^3
    r13_5 = r13_3*r13^2
    r23_5 = r23_3*r23^2

    # Avoid calculating multiple times...
    omm = 1.0 - μ

    # Due to symmetry of second order partials, only 6 (possibly) unique values to calculate

    # Diagonals (share same equation structure)
    Ωxx = 1 - omm / r13_3 - μ / r23_3 + 3omm * (pv[1] + μ)^2 / r13_5 + 3μ * (pv[1] - omm)^2 / r23_5
    Ωyy = 1 - omm / r13_3 - μ / r23_3 + 3omm * (pv[2])^2 / r13_5 + 3μ * (pv[2])^2 / r23_5
    Ωzz =   - omm / r13_3 - μ / r23_3 + 3omm * (pv[3])^2 / r13_5 + 3μ * (pv[3])^2 / r23_5

    # Off-diagonals (share same equation structure)
    Ωxy = 3omm * (pv[1] + μ) * pv[2] / r13_5 + 3μ * (pv[1] - omm) * pv[2] / r23_5
    Ωxz = 3omm * (pv[1] + μ) * pv[3] / r13_5 + 3μ * (pv[1] - omm) * pv[3] / r23_5
    Ωyz = 3omm * (pv[2]) * pv[3] / r13_5 + 3μ * (pv[2]) * pv[3] / r23_5

    # Give back the entire 6x6 matrix even though only first quadrant has values.
    # This is done to avoid suprise as the gradient is a vector of length 6 and the least suprise would
    # be partials with respect to all 6 state variables.
    # A sparse matrix optimization *could* be applied here, but I doubt it would alter much.
    # It's 2020 so I think that 75% sparsity on a 216 byte matrix is a wash...
    # P.S.: Official partial matrix of the 4th of July...?
    @SMatrix [
        Ωxx Ωxy Ωxz 0.0 0.0 0.0;
        Ωxy Ωyy Ωyz 0.0 0.0 0.0;
        Ωxz Ωyz Ωzz 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0
    ]
end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                        JACOBI CONSTANT                                           #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    jacobi_constant(::CrtbpSystem, ::PositionVelocity)

Calculate the Jacobi constant of a state in the CRTBP.

See also: [`jacobi_constant_gradient`](@ref)
"""
function jacobi_constant(sys::CrtbpSystem, pv::PositionVelocity)
    v = pv_velocity_mag(pv)
    Ω = pseudopotential(sys, pv)
    2Ω-v^2
end
"""
    jacobi_constant_gradient(::CrtbpSystem, ::PositionVelocity)

Calculate the gradient of Jacobi constant with respect to all 6 position-velocity variables.

See also: [`jacobi_constant`](@ref)
"""
function jacobi_constant_gradient(sys::CrtbpSystem, pv::PositionVelocity)
    Ωx = pseudopotential_gradient(sys, pv)
    v_vec = pv_velocity(pv)
    vx = @SVector [0.0, 0.0, 0.0, v_vec[1], v_vec[2], v_vec[3]]
    2Ωx .- 2 .* vx
end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                       EQUATIONS OF MOTION                                        #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    crtbp_eom(pv::AbstractArray, sys::CrtbpSystem, t)

Evaluate the equations of motion for the Circular Restricted Three Body Problem (CRTBP) at the
state specified by `pv` in the system specified by `sys`.
The time input, `t`, while required due to compatibility with `DifferentialEquations.jl`, is not used.

Note, here, `pv` can be any `AbstractArray` due to compatibility with `DifferentialEquations`.
However, if a `PositionVelocity` is **not** constructible via `PositionVelocity(pv)`, an error will
be thrown.

See also: [`crtbp_jacobian`](@ref), [`PositionVelocity`](@ref)
"""
function crtbp_eom(q::AbstractArray, sys::CrtbpSystem, t)
    μ = mass_ratio(sys)
    pv = PositionVelocity(q)
    vel   = pv_velocity(pv)
    Ωx    = pseudopotential_gradient(sys, pv)
    @SVector [
        vel[1],
        vel[2],
        vel[3],
        2vel[2] + Ωx[1],
        -2vel[1] + Ωx[2],
        Ωx[3],
    ]
end

"""
    crtbp_jacobian(pv::AbstractArray, sys::CrtbpSystem, t)

Evaluate the jacobian of equations of motion for the Circular Restricted Three Body Problem (CRTBP)
at the state specified by `pv` in the system specified by `sys`.
The time input, `t`, while required due to compatibility with `DifferentialEquations.jl`, is not used.

Note, here, `pv` can be any `AbstractArray` due to compatibility with `DifferentialEquations`.
However, if a `PositionVelocity` is **not** constructible via `PositionVelocity(pv)`, an error will
be thrown.

See also: [`crtbp_eom`](@ref), [`PositionVelocity`](@ref)
"""
function crtbp_jacobian(q::AbstractArray, sys::CrtbpSystem, t)
    pv = PositionVelocity(q)
    Ωij = pseudopotential_hessian(sys, pv)

    Ωxx = Ωij[1,1]
    Ωxy = Ωij[1,2]
    Ωxz = Ωij[1,3]
    Ωyy = Ωij[2,2]
    Ωyz = Ωij[2,3]
    Ωzz = Ωij[3,3]

    @SMatrix [
        0.0  0.0  0.0  1.0  0.0  0.0;
        0.0  0.0  0.0  0.0  1.0  0.0;
        0.0  0.0  0.0  0.0  0.0  1.0;
        Ωxx  Ωxy  Ωxz  0.0  2.0  0.0;
        Ωxy  Ωyy  Ωyz -2.0  0.0  0.0;
        Ωxz  Ωyz  Ωzz  0.0  0.0  0.0
    ]
end
