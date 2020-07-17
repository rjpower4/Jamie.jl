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
struct CrtbpP1{T <: AbstractBody} <: CrtbpPrimary
    body::T
end

"""
    CrtbpP1(body=NullBody())

Construct a new `CrtbpP1` instance.
The internal `body` field defaults to the `NullBody`.

See also: [`CelestialBody`](@ref), [`NullBody`](@ref), ['CrtbpP2'](@ref)
"""
function CrtbpP1(body=NullBody())
    CrtbpP1(body)
end


"""
    CrtbpP2

Empty struct identifying the second (usually smaller) primary in the CRTBP.

See also: [`CrtbpP1`](@ref)
"""
struct CrtbpP2{T <: AbstractBody} <: CrtbpPrimary
    body::T
end

"""
    CrtbpP2(body=NullBody())

Construct a new `CrtbpP2` instance.
The internal `body` field defaults to the `NullBody`.

See also: [`CelestialBody`](@ref), [`NullBody`](@ref), ['CrtbpP1'](@ref)
"""
function CrtbpP2(body=NullBody())
    CrtbpP2(body)
end

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                     CRTBP SYSTEM DEFINITIONS                                     #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    CrtbpSystem{T, CM, CL, CT}(μ::T, name::String, char_mass::CM, char_length::CL, char_time::CT)

Struct defining the parameters for the Crtbp System. Contains the mass ratio parameter, ``\\mu``,
as well as the system name and characteristic quantities.

See also: [`mass_ratio`](@ref), [`characteristic_mass`](@ref), [`characteristic_time`](@ref),
[`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
struct CrtbpSystem{T, D}
    μ::T
    name::String
    dimset::D

    function CrtbpSystem(μ::T, name::String, dimset::D) where {T, D}
        if μ <= 0.0
            throw(DomainError(μ, "Non-positive μ value"))
        end
        new{T, D}(μ, name, dimset)
    end
end

"""
    name_string(::CrtbpSystem)

Retrieve the name of the CRTBP system.

See also:
[`CrtbpSystem`](@ref),
"""
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
function CrtbpSystem(μ::T; name="UNNAMED SYSTEM", dimset::D=nothing) where {T, D}
    CrtbpSystem(μ, name, dimset)
end


"""
    mass_ratio(::CrtbpSystem)

Retrieve the mass ratio, μ, of the specified system.

See also: [`CrtbpSystem`](@ref), [`characteristic_mass`](@ref), [`characteristic_time`](@ref),
[`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
mass_ratio(s::CrtbpSystem) = s.μ

"""
    dimensional_set(::CrtbpSystem{T, D <: DimensionalSet})

Retrieve the dimensional set for the specified system.
"""
dimensional_set(s::CrtbpSystem{T, D}) where {T, D <: DimensionalSet} = s.dimset

"""
    characteristic_mass(::CrtbpSystem)

Retrieve the characteristic mass of the CRTBP system specified. This mass is usually defined as
the sum of the masses of the two primaries.

See also: [`CrtbpSystem`](@ref), 
[`mass_ratio`](@ref), 
[`dimensional_set`](@ref),
[`characteristic_time`](@ref),
[`characteristic_length`](@ref), 
[`characteristic_velocity`](@ref),
[`characteristic_acceleration`](@ref)
"""
characteristic_mass(s::CrtbpSystem) = dimensional_set(s) |> characteristic_mass

"""
    characteristic_length(::CrtbpSystem)

Retrieve the characteristic length of the CRTBP system specified. This length is usually defined as
the semi-major axis of the circular orbit of P2 in the CRTBP.

See also: [`CrtbpSystem`](@ref), 
[`dimensional_set`](@ref),
[`mass_ratio`](@ref), [`characteristic_mass`](@ref),
[`characteristic_time`](@ref), 
[`characteristic_acceleration`](@ref)
"""
characteristic_length(s) = dimensional_set(s) |> characteristic_length

"""
    characteristic_length(::CrtbpSystem)

Retrieve the characteristic time of the CRTBP system specified. This time is usually defined as the
time such that the universal gravitational constant in the non-dimensional system is unity.

See also: [`CrtbpSystem`](@ref), 
[`dimensional_set`](@ref),
[`mass_ratio`](@ref), [`characteristic_mass`](@ref),
[`characteristic_length`](@ref), 
[`characteristic_velocity`](@ref),
[`characteristic_acceleration`](@ref)
"""
characteristic_time(s::CrtbpSystem) = dimensional_set(s) |> characteristic_time

"""
    characteristic_velocity(::CrtbpSystem)

Retrieve the characteristic velocityof the CRTBP system specified. This is equivalent to the
characteristic length divided by the characteristic time.

See also: [`CrtbpSystem`](@ref), 
[`dimensional_set`](@ref),
[`mass_ratio`](@ref), [`characteristic_mass`](@ref),
[`characteristic_length`](@ref),
[`characteristic_acceleration`](@ref)
"""
characteristic_velocity(s::CrtbpSystem) = dimensional_set(s) |> characteristic_velocity


"""
    characteristic_acceleration(::CrtbpSystem)

Retrieve the characteristic accerlation of the CRTBP system specified. This is equivalent to the
characteristic length divided by the characteristic time squared.

See also: [`CrtbpSystem`](@ref), 
[`dimensional_set`](@ref),
[`mass_ratio`](@ref), [`characteristic_mass`](@ref),
[`characteristic_length`](@ref),
[`characteristic_velocity`](@ref)
"""
characteristic_acceleration(s::CrtbpSystem) = dimensional_set(s) |> characteristic_acceleration

"""
    dimensionalize(::CrtbpSystem{T,T}, ::PositionVelocity)

Dimensionalize the position and velocity components using the characteristic quantities associated
with the CRTBP system.

See also: [`nondimensionalize`](@ref)
[`CrtbpSystem`](@ref), [`characteristic_length`](@ref), [`characteristic_velocity`](@ref)
"""
function dimensionalize(sys::CrtbpSystem, pv::PositionVelocity)
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
function nondimensionalize(sys::CrtbpSystem, pv::PositionVelocity)
    cl = characteristic_length(sys)
    cv = characteristic_velocity(sys)
    scale_vec = @SVector [cl, cl, cl, cv, cv, cv]
    pv ./ scale_vec
end

"""
    PositionVelocity(::CrtbpSystem, ::CrtbpP1)

Return the position and velocity of ``P_1`` in the CRTBP.
"""
function PositionVelocity(sys::CrtbpSystem, ::CrtbpP1)
    PositionVelocity(
        -mass_ratio(sys),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    )
end

"""
    PositionVelocity(::CrtbpSystem, ::CrtbpP2)

Return the position and velocity of ``P_2`` in the CRTBP.
"""
function PositionVelocity(sys::CrtbpSystem, ::CrtbpP2)
    PositionVelocity(
        1.0 - mass_ratio(sys),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    )
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
    CrtbpSystem(pv::AbstractArray, p, t)

Evaluate the equations of motion for the Circular Restricted Three Body Problem (CRTBP) at the
state specified by `pv`.
The time input, `t`, and the parameter input `p`, 
while required due to compatibility with `DifferentialEquations.jl`, are not used.

Note, here, `pv` can be any `AbstractArray` due to compatibility with `DifferentialEquations`.
However, if a `PositionVelocity` is **not** constructible via `PositionVelocity(pv)`, an error will
be thrown.

See also: [`CrtbpSystem`](@ref), [`PositionVelocity`](@ref)
"""
function (sys::CrtbpSystem)(q::AbstractArray, p, t)
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

function (sys::CrtbpSystem)(q::AbstractArray, t)
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

function (sys::CrtbpSystem)(q::AbstractArray)
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


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                  LAGRANGE POINT DETERMINATION                                    #
# ************************************************************************************************ #
# ************************************************************************************************ #
module LagrangeHelper

export __crtbp_eq_init_guesses
export __crtbp_eq_helper1, __crtbp_eq_helper2, __crtbp_eq_helper3

function __crtbp_eq_init_guesses(μ)
    η = (μ / 3.0)^(1.0/3.0)
    σ = 7.0μ/12.0

    xl1 = 1.0 - μ -
        η * (1.0 - η/3.0 - η^2/9.0 - 23.0η^3/81.0 + 151.0η^4/243.0 - η^5/9.0)
    xl2 = 1.0 - μ +
        η * (1.0 + η/3.0 - η^2/9.0 - 31.0η^3/81.0 - 119.0η^4/243.0 - η^5/9.0)
    xl3 = -μ - 1.0 +
        σ * (1.0 + 23.0σ^2/84.0 + 23.0σ^3/84.0 + 761.0σ^4/2352.0 +
             3163.0σ^5/7056.0 + 30703.0σ^6/49392.0)
    (xl1, xl2, xl3)
end



# ************************************************************************************************ #
# ************************************************************************************************ #
#                           HELPER FUNCTIONS FOR NR                            #
# ************************************************************************************************ #
# ************************************************************************************************ #

function __crtbp_eq_helper1(μ, γ)
    ( 1.0 - μ - γ - (1.0 - μ) / (1.0 - γ)^2 + μ / γ^2,
     -1.0 - 2.0 * (1.0 - μ) / (1.0 - γ)^3 - 2.0μ / γ^3)
end

function __crtbp_eq_helper2(μ, γ)
    ((1.0 - μ) / (1.0 + γ)^2 + μ / γ^2 - 1.0 + μ - γ,
    -2.0 * (1.0 - μ) / (1.0 + γ)^3 - 2.0 * μ / γ^3 - 1.0)
end


function __crtbp_eq_helper3(μ, γ)
    (-μ - γ + (1 - μ) / γ^2 + μ / (γ + 1.0)^2,
     -1.0 - 2.0*(1.0 - μ) / γ^3 - 2.0μ / (γ + 1.0)^3)
end


function __crtbp_eq_nr(F, μ, x0, tolerance, max_iter)
    f, df = F(μ, x0)
    for i = 1:max_iter
        if abs(f) < tolerance
            conv = true
            break
        else
            x0 -= f/df
            f, df = F(μ, x0)
        end
    end
    x0
end

end # module LagrangeHelper


"""
    equilibrium_solutions(::CrtbpSystem; [tolerance=1e-12], [max_iter=20])

Obtain the (5) equilibrium solutions for the CRTBP often called Lagrange points.
These are returned as states. The returned values are not guaranteed to have
converged, but with default options the likelihood of non-convergence is near zero.
"""
function equilibrium_solutions(sys::CrtbpSystem; tolerance=1e-12, max_iter=20)

    μ = mass_ratio(sys)

    # Get the initial guesses for the collinear points
    x1_guess, x2_guess, x3_guess = LagrangeHelper.__crtbp_eq_init_guesses(μ)

    # Convert the x coordinate guesses to the transformed γ values
    γ1_guess = 1.0 - μ - x1_guess
    γ2_guess = x2_guess - 1.0 + μ
    γ3_guess = -μ - x3_guess

    # Run newton raphson for the three different
    γ1 = LagrangeHelper.__crtbp_eq_nr(LagrangeHelper.__crtbp_eq_helper1, μ, γ1_guess, tolerance, max_iter)
    γ2 = LagrangeHelper.__crtbp_eq_nr(LagrangeHelper.__crtbp_eq_helper2, μ, γ2_guess, tolerance, max_iter)
    γ3 = LagrangeHelper.__crtbp_eq_nr(LagrangeHelper.__crtbp_eq_helper3, μ, γ3_guess, tolerance, max_iter)

    # Convert the (hopefully) converged γ's to x values. This is the most
    # obvious place to see the transform as well.
    x1 = 1.0 - μ - γ1
    x2 = 1.0 - μ + γ2
    x3 = 0.0 - μ - γ3

    L1 = PositionVelocity(x1, 0.0, 0.0, 0.0, 0.0, 0.0)
    L2 = PositionVelocity(x2, 0.0, 0.0, 0.0, 0.0, 0.0)
    L3 = PositionVelocity(x3, 0.0, 0.0, 0.0, 0.0, 0.0)

    # L4 and L5 have the same analytically determined solution for all systems
    l45_y = sqrt(3.0) / 2.0
    L4 = PositionVelocity(0.5-μ,  l45_y, 0.0, 0.0, 0.0, 0.0)
    L5 = PositionVelocity(0.5-μ, -l45_y, 0.0, 0.0, 0.0, 0.0)

    (L1, L2, L3, L4, L5)
end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                       ROTATION ALGORITHMS                                        #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    shift_basepoint(::CrtbpSystem, ::CrtbpPrimary, ::PositionVelocity)

Return the position velocity transformed so that it is with respect to the basepoint specified by the primary.
Note, this is still in the *rotating frame*.

See also:
[`CrtbpSystem`](@ref),
[`CrtbpP1`](@ref),
[`CrtbpP2`](@ref),
[`PositionVelocity`](@ref),
[`shift_to_inertial`](@ref)
"""
function shift_basepoint(sys::CrtbpSystem, primary::CrtbpPrimary, pv::PositionVelocity)
    pv - PositionVelocity(sys, primary)
end

"""
    shift_to_inertial(::CrtbpSystem, ::CrtbpPrimary, ::PositionVelocity; t, [initial_phase=0])

Translate and rotate the input such that it has a basepoint at the specified primary and is in the 
"inertial" frame. 
This inertial frame is still that of the idealized CRTBP, i.e., it is a simple 2π-periodic affine transformation and
no ephemerides are incorporated.

See also:
[`CrtbpSystem`](@ref),
[`CrtbpP1`](@ref),
[`CrtbpP2`](@ref),
[`PositionVelocity`](@ref),
[`shift_basepoint`](@ref)
"""
function shift_to_inertial(sys::CrtbpSystem, primary::CrtbpPrimary, pv::PositionVelocity; t, initial_phase=0.0)
    ϕ = t  + initial_phase
    C = @SMatrix [
         cos(ϕ) -sin(ϕ) 0.0 0.0     0.0    0.0;
         sin(ϕ)  cos(ϕ) 0.0 0.0     0.0    0.0;
         0.0     0.0    1.0 0.0     0.0    0.0;
        -sin(ϕ) -cos(ϕ) 0.0 cos(ϕ) -sin(ϕ) 0.0;
         cos(ϕ) -sin(ϕ) 0.0 sin(ϕ)  cos(ϕ) 0.0;
         0.0     0.0    0.0 0.0     0.0    1.0
    ]
    C * shift_basepoint(sys, primary, pv)
end