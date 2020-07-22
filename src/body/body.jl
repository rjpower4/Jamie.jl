include("potential.jl")
include("shape.jl")

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                      ABSTRACT BODY DEFINION                                      #
# ************************************************************************************************ #
# ************************************************************************************************ #

"""
    AbstractBody

Abstract base type for all "body-like" objects. Bodies in Jamie are amorphous in design. The most
common implementation is the [`CelestialBody`](@ref) which defines a planet/moon/asteroid type
object. The method defined on the `AbstractBody` is the [`name_string`](@ref) method which returns
"UNNAMED BODY" if not redifined for a specific type.

See also:
[`CelestialBody`](@ref),
[`NullCelestialBody`](@ref)
"""
abstract type AbstractBody end

"""
    name_string(::AbstractBody)

Return the name of the body as a string.
"""
name_string(::AbstractBody) = "UNNAMED BODY"


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                 ABSTRACT CELESTIAL BODY STRUCTURE                                #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    AbstractCelestialBody <: AbstractBody

This subtype of [`AbstractBody`](@ref) specifies that the body is a "Celestial" body like a
planet, moon, or asteroid.
"""
abstract type AbstractCelestialBody <: AbstractBody end

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                   NULL CELESTIAL BODY STRUCTURE                                  #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    NullCelestialBody <: AbstractBody

This is an empty celestial body type.
"""
struct NullCelestialBody <: AbstractCelestialBody end

"""
    name_string(::NullCelestialBody)

Returns the name of the null body, always "NULL BODY".
"""
name_string(::NullCelestialBody) = "NULL BODY"

"""
    shape_model(::NullCelestialBody)

Returns a null shape model.

See also:
[`NullShapeModel`](@ref)
"""
shape_model(::NullCelestialBody) = NullShapeModel()

"""
    potential_model(::NullCelestialBody)

Returns a null potential model.

See also:
[`NullPotentialModel`](@ref)
"""
potential_model(::NullCelestialBody) = NullPotential()

# ************************************************************************************************ #
# ************************************************************************************************ #
#                                      CELESTIAL BODY STRUCTURE                                    #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    CelestialBody{T, I <: Integer}(name::String, gm::T, spice_id::I, mean_radius::T) <: AbstractCelestialBody

Struct defining a celestial body (e.g. Moon, Planet, Asteroid).

See also:
[`PointMassPotential`](@ref),
[`SphericalShapeModel`](@ref),
[`NullCelestialBody`](@ref)
"""
struct CelestialBody{I <: Integer, P <: PointMassPotential, S <: AbstractShapeModel} <: AbstractCelestialBody
    name::String
    spice_id::I
    potential::P
    shape::S
end

"""
    name_string(cb::CelestialBody)

Retrieve the name of the celestial body.

See also:
[`CelestialBody`](@ref),
[`spice_identifier`](@ref),
[`potential_model`](@ref),
[`shape_model`](@ref),
[`mean_radius`](@ref),
[`gravitational_parameter`](@ref)
"""
name_string(cb::CelestialBody) = cb.name

"""
    spice_identifier(cb::CelestialBody)

Retrieve the integer SPICE identifier of the celestial body.
More information on SPICE ID codes can be found
[here](http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/06_naif_ids.pdf).

See also:
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`potential_model`](@ref),
[`shape_model`](@ref),
[`mean_radius`](@ref),
[`gravitational_parameter`](@ref)
"""
spice_identifier(cb::CelestialBody) = cb.spice_id

"""
    potential_model(::CelestialBody)

Retrieve the potential modlel of the body.

See also:
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`spice_identifier`](@ref),
[`shape_model`](@ref),
[`mean_radius`](@ref),
[`gravitational_parameter`](@ref)
"""
potential_model(cb::CelestialBody) = cb.potential

"""
    shape_model(::CelestialBody)

Retrieve the shape model of the body.

See also:
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`spice_identifier`](@ref),
[`potential_model`](@ref),
[`mean_radius`](@ref),
[`gravitational_parameter`](@ref)
"""
shape_model(cb::CelestialBody) = cb.shape

"""
    mean_radius(::CelestialBody)

Retrieve the mean radius of the body if a non-null shape model is defined.
The particular calculation of the mean radius is implemented by the shape model.

See also:
[`SphericalShapeModel`](@ref),
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`spice_identifier`](@ref),
[`potential_model`](@ref),
[`shape_model`](@ref),
[`gravitational_parameter`](@ref)
"""
mean_radius(cb::CelestialBody) = shape_model(cb) |> mean_radius

"""
    gravitational_parameter(::CelestialBody)

Retrieve the gravitational parameter of the body if a non-null potential model is defined.
The particular calculation of the gravitational parameter is implemented by the potential model.

See also:
[`PointMassPotential`](@ref),
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`spice_identifier`](@ref),
[`potential_model`](@ref),
[`shape_model`](@ref),
[`mean_radius`](@ref),
"""
gravitational_parameter(cb::CelestialBody) = potential_model(cb) |> gravitational_parameter

"""
    gravitational_potential(::CelestialBody)

Retrieve the gravitational potential of the body if a non-null potential model is defined.
The particular calculation of the gravitational potential is implemented by the potential model.

See also:
[`PointMassPotential`](@ref),
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`spice_identifier`](@ref),
[`potential_model`](@ref),
[`shape_model`](@ref),
[`mean_radius`](@ref),
"""
gravitational_potential(cb::CelestialBody, q) = gravitational_potential(potential_model(cb), q)

"""
    gravitational_acceleration(::CelestialBody, r)

Retrieve the gravitational acceleration of the body if a non-null potential model is defined.
The value, `r`, can be a scalar distance as well as a vector (`Array`, `PositionVelocity`, ...).
If a scalar is passed, the acceleration magnitude is returned.

See also:
[`PointMassPotential`](@ref),
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`spice_identifier`](@ref),
[`potential_model`](@ref),
[`shape_model`](@ref),
[`mean_radius`](@ref),
"""
gravitational_acceleration(cb::CelestialBody, q) = gravitational_acceleration(potential_model(cb), q)
