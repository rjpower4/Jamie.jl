# ************************************************************************************************ #
# ************************************************************************************************ #
#                                     ABSTRACT POTENTIAL MODEL                                     #
# ************************************************************************************************ #
# ************************************************************************************************ #
abstract type AbstractPotential end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                       NULL POTENTIAL MODEL                                       #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    NullPotential

Empty potential model used as an indicator that no potential is defined.

See also: [`PointMassPotential`](@ref)
"""
struct NullPotential <: AbstractPotential end


# ************************************************************************************************ #
# ************************************************************************************************ #
#                                    POINT MASS POTENTIAL MODEL                                    #
# ************************************************************************************************ #
# ************************************************************************************************ #
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

"""
    gravitational_potential(::PointMassPotential, r)

Calculate the gravitational potential of another massive body at
radius, `r`, with respect to the center of mass of the potential model.
"""
gravitational_potential(pmp::PointMassPotential, r) = -gravitational_parameter(pmp) / r

"""
    gravitational_potential(::PointMassPotential, r::AbstractArray)

Calculate the gravitational potential of another massive body with
radius vector, `r`, with respect to the center of mass of the potential model.
"""
function gravitational_potential(pmp::PointMassPotential, r::AbstractArray)
    gravitational_potential(pmp, sqrt(r[1]^2 + r[2]^2 + r[3]^2))
end

"""
    gravitational_potential(::PointMassPotential, pv::PositionVelocity)

Calculate the gravitational potential of another massive body with
state, `pv`, with respect to the center of mass of the potential model.
"""
function gravitational_potential(pmp::PointMassPotential, pv::PositionVelocity)
    gravitational_potential(pmp, position_mag(pv))
end

"""
    gravitational_acceleration(::PointMassPotential, r)

Calculate the acceleration magnitude due to gravity for a body at position radius, `r`,
with respect to the point mass.
"""
function gravitational_acceleration(pmp::PointMassPotential, r)
    gravitational_parameter(pmp) / r^2
end

"""
    gravitational_acceleration(::PointMassPotential, r::AbstractArray)

Calculate the acceleration vector due to gravity for a body at position vector, `r`,
with respect to the point mass.
"""
function gravitational_acceleration(pmp::PointMassPotential, r::AbstractArray)
    rmag = norm(r)
    accel_mag = gravitational_acceleration(pmp, rmag)
    -accel_mag / rmag * r
end

"""
    gravitational_acceleration(::PointMassPotential, pv::PositionVelocity)

Calculate the acceleration vector due to gravity for a body with position and velocity
with respect to the point mass defined by `pv`.
"""
function gravitational_acceleration(pmp::PointMassPotential, pv::PositionVelocity)
    gravitational_acceleration(pmp, position(pv))
end
