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
abstract type CelestialBody <: AbstractCelestialBody end

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
name_string(cb::CelestialBody) = "Celestial Body"

"""
    spice_identifier(cb::CelestialBody)

Retrieve the integer SPICE identifier of the celestial body.
More information on SPICE ID codes can be found
[here](http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/06_naif_ids.pdf).

See also:
[`CelestialBody`](@ref),
[`name_string`](@ref),
[`mean_radius`](@ref),
[`gravitational_parameter`](@ref)
"""
spice_identifier(cb::CelestialBody) = throw(MethodError(spice_identifier, cb))

"""
    mean_radius(::CelestialBody)

Retrieve the mean radius of the body if a non-null shape model is defined.
The particular calculation of the mean radius is implemented by the shape model.
"""
mean_radius(cb::CelestialBody) = throw(MethodError(mean_radius, cb))

equatorial_radius(cb::CelestialBody) = throw(MethodError(equatorial_radius, cb))

polar_radius(cb::CelestialBody) = throw(MethodError(polar_radius, cb))

ellipticity(cb::CelestialBody) = throw(MethodError(ellipticity, cb))

mean_density(cb::CelestialBody) = throw(MethodError(mean_density, cb))

number_natural_satellites(cb::CelestialBody) = throw(MethodError(number_natural_satellites, cb))

mean_semimajor_axis(cb::CelestialBody) = throw(MethodError(mean_semimajor_axis, cb))

mean_eccentricity(cb::CelestialBody) = throw(MethodError(mean_eccentricity, cb))

solar_irradiance(cb::CelestialBody) = throw(MethodError(solar_irradiance, cb))

parent_body(cb::CelestialBody) = throw(MethodError(parent_body, cb))

"""
    gravitational_parameter(::CelestialBody)

Retrieve the gravitational parameter of the body if a non-null potential model is defined.
The particular calculation of the gravitational parameter is implemented by the potential model.
"""
gravitational_parameter(cb::CelestialBody) = throw(MethodError(gravitational_parameter, cb))

total_mass(cb::CelestialBody) = throw(MethodError(total_mass, cb))
total_volume(cb::CelestialBody) = throw(MethodError(total_volume, cb))