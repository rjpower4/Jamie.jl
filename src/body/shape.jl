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
