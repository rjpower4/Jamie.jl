# ************************************************************************************************ #
# ************************************************************************************************ #
#                                        POTENTIAL MODELS                                          #
# ************************************************************************************************ #
# ************************************************************************************************ #
abstract type AbstractPotential end

"""
    NullPotential

Empty potential model used as an indicator that no potential is defined.

See also: [`PointMassPotential`](@ref)
"""
struct NullPotential <: AbstractPotential end

"""
    PointMassPotential(gm::T)

Simple potential model for a point mass (or centrobaric) body of gravitational parameter `gm`.

See also: [`NullPotential`](@ref)
"""
struct PointMassPotential{T}
    gm::T
end

"""
    gravitational_parameter(::PointMassPotential)

Retrieve the gravitational parameter, `GM`, of the potential model.
This is the standard 2-body μ value.

See also: [`PointMassPotential`](@ref)
"""
gravitational_parameter(pmp::PointMassPotential) = pmp.gm


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                           SHAPE MODELS                                           #
# ************************************************************************************************ #
# ************************************************************************************************ #
abstract type AbstractShapeModel end

"""
    NullShapeModel
Empty shape model used as an indicator that no shape is defined.

See also: [`SphericalShapeModel`](@ref)
"""
struct NullShapeModel <: AbstractShapeModel end

"""
    SphericalShapeModel(radius::T)

Simple shape model describing a sphere of radius, `radius`.

See also: [`NullShapeModel`](@ref)
"""
struct SphericalShapeModel{T} <: AbstractShapeModel
    radius::T
end

"""
    mean_radius(::SphericalShapeModel)

Retrieve the mean radius of the spherical shape model.
Clearly for this model, the mean radius is simple equal to the radius, `radius`.

See also: [`NullShapeModel`](@ref)
"""
mean_radius(m::SphericalShapeModel) = m.radius


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                      CELESTIAL BODY STRUCTURE                                    #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    CelestialBody{T, I <: Integer}(name::String, gm::T, spice_id::I, mean_radius::T)

Struct defining a celestial body (e.g. Moon, Planet, Asteroid).

See also: [`PointMassPotential`](@ref), [`SphericalShapeModel`](@ref)
"""
struct CelestialBody{I <: Integer, P <: PointMassPotential, S <: AbstractShapeModel}
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
