"""
    DEFAULT_INTEGRATION_ALGORITHM

The default integration method used for `DifferentialEquations.jl` integrations. A full list of
methods can be found [here](https://diffeq.sciml.ai/stable/solvers/ode_solve/)
"""
const DEFAULT_INTEGRATION_ALGORITHM = Vern9

"""
    DEFAULT_SOLVER_OPTIONS

The default solver options to be passed to the `DifferentialEquations.jl` solver. Currently, only
the absolute and relative tolerances are included. A full list of solver options can be found
[here](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).
"""
const DEFAULT_SOLVER_OPTIONS = Dict(
    :abstol => 1e-12, # Default Absolute Tolerance
    :reltol => 1e-12, # Default Relative Tolerance
)
