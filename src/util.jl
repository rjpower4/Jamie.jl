# ************************************************************************************************ #
# ************************************************************************************************ #
#                                        POINTS ON SPHERE                                          #
# ************************************************************************************************ #
# ************************************************************************************************ #
"""
    points_on_sphere(n; [radius=1.0], [center=(0.0, 0.0, 0.0)])

Generate a set of points on a sphere *approximately* uniformly distributed via the Fibonacci lattice
method.
The points are returned as an `Array` of `n` `SVector`s of length 3.
A specific `radius` and `center` can be optionally given as key-word arguments.
"""
function points_on_sphere(n; radius=1.0, center=(0.0, 0.0, 0.0))
    output = Array{SVector{3, Float64}, 1}(undef, n)
    gr = (1.0 + sqrt(5)) ./ 2.0
    for k in 1:n
        i = (k - 1.0) + 0.5
        ϕ = acos(1.0 - 2i / n)
        θ = 2π * i  / gr
        output[k] = SVector(
            radius * cos(θ) * sin(ϕ) + center[1],
            radius * sin(θ) * sin(ϕ) + center[2],
            radius * cos(ϕ) + center[3]
        )
    end
    output
end
