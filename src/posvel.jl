"""
    PositionVelocity{T} <: FieldVector{6, T}

Struct representing the position-velocity vector.
`T` designates the types of component in the vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
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
    position(::PositionVelocity)

Get the position components of the position-velocity vector.

See also: [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
Base.position(pv::PositionVelocity) = @SVector [ pv.x, pv.y, pv.z ]

"""
    position_magnitude(::PositionVelocity)

Return the magnitude of the position components of the position-velocity vector.
This is equivalent to the distance from the origin of the position vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
position_magnitude(pv::PositionVelocity) = (norm ∘ position)(pv)

"""
    position_direction(::PositionVelocity)

Return the unit vector parallel to the position sub-vector in the specified position-velocity vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
position_direction(pv::PositionVelocity) = position(pv) ./ position_magnitude(pv)

"""
    velocity(::PositionVelocity)

Get the velocity components of the position-velocity vector.

See also: [`position`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
[`velocity_direction`](@ref)
"""
velocity(pv::PositionVelocity) = @SVector [ pv.vx, pv.vy, pv.vz ]

"""
    velocity_magnitude(::PositionVelocity)

Return the magnitude of the velocity components of the position-velocity vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref),
[`velocity_direction`](@ref)
"""
velocity_magnitude(pv::PositionVelocity) = (norm ∘ velocity)(pv)

"""
    velocity_direction(::PositionVelocity)

Return the unit vector parallel to the velocity sub-vector in the specified position-velocity vector.

See also: [`position`](@ref), [`velocity`](@ref),
[`position_magnitude`](@ref), [`position_direction`](@ref), [`velocity_magnitude`](@ref),
"""
velocity_direction(pv::PositionVelocity) = velocity(pv) ./ velocity_magnitude(pv)


"""
    dualize(pv::PositionVelocity)

Build a dual position-velocity to enable evaluation of partial derivatives.
"""
function ForwardDiff.dualize(pv::PositionVelocity{T}) where {T}
    PositionVelocity(
        ForwardDiff.dualize(PositionVelocity{T}, pv)
    )
end
