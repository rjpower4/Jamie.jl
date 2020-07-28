abstract type AbstractTwoBodyElementSet{N, T} <: FieldVector{N, T} end

"""
    KeplerianElementSet{T}

Structure containing the six standard Keplerian elements .
Uses true anomaly, ``\\theta``, as the "fast variable", so any other
anomalies requested are calculated (read: overhead).

See also:
[`semimajor_axis`](@ref),
[`eccentricity`](@ref),
[`inclination`](@ref),
[`argument_of_periapsis`](@ref),
[`right_ascension`](@ref),
[`true_anomaly`](@ref)
"""
struct KeplerianElementSet{T} <: AbstractTwoBodyElementSet{6, T}
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    θ::T
end

"""
    semimajor_axis(::KeplerianElementSet)

Retrieve the semimajor axis of the orbit specified.

See also:
[`KeplerianElementSet`](@ref),
[`eccentricity`](@ref),
[`inclination`](@ref),
[`argument_of_periapsis`](@ref),
[`right_ascension`](@ref),
[`true_anomaly`](@ref)
"""
semimajor_axis(ke::KeplerianElementSet) = ke.a

"""
    eccentricity(::KeplerianElementSet)

Retrieve the eccentricity of the orbit specified.

See also:
[`KeplerianElementSet`](@ref),
[`semimajor_axis`](@ref),
[`eccentricity`](@ref),
[`inclination`](@ref),
[`argument_of_periapsis`](@ref),
[`right_ascension`](@ref),
[`true_anomaly`](@ref)
"""
eccentricity(ke::KeplerianElementSet) = ke.e

"""
    inclination(::KeplerianElementSet)

Retrieve the inclination of the orbit specified.

See also:
[`KeplerianElementSet`](@ref),
[`semimajor_axis`](@ref),
[`eccentricity`](@ref),
[`argument_of_periapsis`](@ref),
[`right_ascension`](@ref),
[`true_anomaly`](@ref)
"""
inclination(ke::KeplerianElementSet) = ke.i

"""
    arguement_of_periapsis(::KeplerianElementSet)

Retrieve the argument of periapsis of the orbit specified.

See also:
[`KeplerianElementSet`](@ref),
[`semimajor_axis`](@ref),
[`eccentricity`](@ref),
[`inclination`](@ref),
[`right_ascension`](@ref),
[`true_anomaly`](@ref)
"""
argument_of_periapsis(ke::KeplerianElementSet) = ke.ω

"""
    right_ascension(::KeplerianElementSet)

Retrieve the right ascension of the orbit specified.

See also:
[`KeplerianElementSet`](@ref),
[`semimajor_axis`](@ref),
[`eccentricity`](@ref),
[`inclination`](@ref),
[`argument_of_periapsis`](@ref),
[`true_anomaly`](@ref)
"""
right_ascension(ke::KeplerianElementSet) = ke.Ω

"""
    true_anomaly(::KeplerianElementSet)

Retrieve the true anomaly of the orbit specified.

See also:
[`KeplerianElementSet`](@ref),
[`semimajor_axis`](@ref),
[`eccentricity`](@ref),
[`inclination`](@ref),
[`argument_of_periapsis`](@ref),
[`right_ascension`](@ref),
"""
true_anomaly(ke::KeplerianElementSet) = ke.θ

"""
    parameter(ke::KeplerianElementSet)

Calculate the parameter (semi-latus rectum) of the Keplerian orbit elements.
"""
parameter(ke::KeplerianElementSet) = semimajor_axis(ke) * (1.0 - eccentricity(ke)^2)

"""
    distance_from_primary(p, e, θ)

Calculate the distance from the primary for orbit state define by parameter, p,
eccentricity, e, and true anomaly, ``\\theta``.
"""
distance_from_primary(p, e, θ) = p / (1.0 + e * cos(θ))

"""
    distance_from_primary(::KeplerianElementSet)

Calculate the distance from the primary for orbit state define by Keplerian elements.
"""
distance_from_primary(ke::KeplerianElementSet) =
    distance_from_primary(parameter(ke), eccentricity(ke), true_anomaly(ke))

"""
    periapsis_distance(ke::KeplerianElementSet)

Calculate the distance at periapsis of the specified orbit.

See also: [`apoapsis_distance`](@ref)
"""
periapsis_distance(ke::KeplerianElementSet) = parameter(ke) * (1.0 - eccentricity(ke))

"""
    apoapsis_distance(ke::KeplerianElementSet)

Calculate the distance at apoapsis of the specified orbit.

See also: [`periapsis_distance`](@ref)
"""
apoapsis_distance(ke::KeplerianElementSet) = parameter(ke) * (1.0 + eccentricity(ke))

"""
    angular_momentum(gm, ke::KeplerianElementSet)

Calculate the angular momentum of the orbit defined by the Keplerian element set around a primary
with the graviational parameter, `gm`.
"""
angular_momentum(gm, ke::KeplerianElementSet) = sqrt(gm * parameter(ke))

"""
    angular_momentum(::AbstractBody, ::KeplerianElementSet)

Calculate the angular momentum of the orbit defined by the Keplerian element set around a primary
define by the `AbstractBody`.
"""
angular_momentum(ab::AbstractBody, ke::KeplerianElementSet) =
    angular_momentum(gravitational_parameter(ab), ke)

"""
    KeplerianElementSet(gm, pv::PositionVelocity)

Convert position-velocity values to Keplerian elements given the gravitational parameter.
The position and velocity components are assumed to be defined in an inertial frame with its 
basepoint at the primary with parameter, `gm`. 
The resultant Keplerian element set will be, therefore, defined with respect to that same frame.
"""
function KeplerianElementSet(gm, pv::PositionVelocity)
    R = position(pv)
    V = velocity(pv)
    r = norm(R)
    v = norm(V)

    spec_energy = 0.5 * v^2 - gm / r
    a = -gm ./ (2spec_energy)

    H = cross(R, V)
    h = norm(H)

    i = acos(H[3] / h)

    eK = @SVector [0.0, 0.0, 1.0]
    nodal_axis = cross(eK, H) ./ h
    nodal_axis_unit = nodal_axis ./ norm(nodal_axis)

    e_vec = (1.0 / gm) .* cross(V, H) .- R ./ r
    e = norm(e_vec)

    if H[2] == 0.0
        nodal_axis_unit = @SVector [cos(0.0), sin(0.0), 0.0]
    end
    Ω = acos(nodal_axis_unit[1])
    if nodal_axis_unit[2] < 0.0
        Ω = 2π - Ω
    end

    ω = acos(dot(nodal_axis_unit, e_vec) ./ e ./ norm(nodal_axis_unit))
    if i != 0.0 && e_vec[3] < 0.0
        ω = 2π - ω
    end
    if i == 0.0 && e_vec[2] < 0.0
        ω = -ω
    end

    θ = 0.0
    if dot(V, R) > 0
        θ = acos(dot(e_vec, R) ./ (r * e))
    else
        θ = 2π - acos(dot(e_vec, R) ./ (r * e))
    end

    KeplerianElementSet(
        a,
        e,
        i,
        ω,
        Ω,
        θ
    )
end

"""
    KeplerianElementSet(::AbstractBody, ::PositionVelocity)

Convert position-velocity values to Keplerian elements given the defined body.
The body is assumed to be centrobaric, i.e., the gm value is simply queried then used directly.
The position and velocity components are assumed to be defined in an inertial frame with its 
basepoint at the center of the body.
The resultant Keplerian element set will be, therefore, defined with respect to that same frame
and basepoint.

See also:
[`CelestialBody`](@ref),
[`gravitational_parameter`](@ref),
"""
function KeplerianElementSet(ab::AbstractBody, pv::PositionVelocity)
    KeplerianElementSet(gravitational_parameter(ab), pv)
end