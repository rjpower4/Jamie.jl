"""
    FlowMap{SYS, SO, P, A}

A `FlowMap` is a structure enabling the flow mapping, or integration, of initial conditions within 
a specified dynamical system. The map contains the dynamical system as well as the default solver 
options, default parameters, and integration algorithm for the integration performed using 
`DifferentialEquations.jl`.

See also:
[`dynamical_system`](@ref),
[`default_solver_options`](@ref),
[`default_parameters`](@ref),
[`solver_algorithm`](@ref),
"""
struct FlowMap{SYS, SO, P, A}
    system::SYS   # Dynamical System
    sopts::SO     # Default Solver Options
    params::P     # Default Parameters
    algorithm::A  # Solver algorithm for DifferentialEquations.jl
end

"""
    dynamical_system(::FlowMap)

Return the dynamical system in use by the flow map.

See also:
[`FlowMap`](@ref),
[`default_solver_options`](@ref),
[`default_parameters`](@ref),
[`solver_algorithm`](@ref),
"""
dynamical_system(fm::FlowMap) = fm.system

"""
    default_solver_options(::FlowMap)


Return the default solver options in use by the flow map.
To see a full list of solver options, see 
[here](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).
These options can be overwritten in the map call.

See also:
[`FlowMap`](@ref),
[`default_solver_options`](@ref),
[`default_parameters`](@ref),
[`solver_algorithm`](@ref),
"""
default_solver_options(fm::FlowMap) = fm.sopts

"""
    default_parameters(::FlowMap)


Return the default parameters in use by the flow map.
These parameters can be overwritten in the map call.

See also:
[`FlowMap`](@ref),
[`default_solver_options`](@ref),
[`default_solver_options`](@ref),
[`solver_algorithm`](@ref),
"""
default_parameters(fm::FlowMap) = fm.params

"""
    solver_algorithm(fm::FlowMap)

Return the solver algorithm used by the flow map.
For a full list of algorithms see [here](https://diffeq.sciml.ai/stable/solvers/ode_solve/).

See also:
[`FlowMap`](@ref),
[`default_solver_options`](@ref),
[`default_solver_options`](@ref),
[`default_parameters`](@ref),
"""
solver_algorithm(fm::FlowMap) = fm.algorithm

"""
    FlowMap(sys; [solve_opts, default_params, algorithm])

Construct a new `FlowMap` with optional arguments initialized with default values.
The default solver options can be seen in `DEFAULT_SOLVER_OPTIONS`.
The default integration algorithm can be seen in 
`DEFAULT_INTEGRATION_ALGORITHM`.
The default parameters are `nothing`, so if a dynamical model requires a specific "null" 
parameter value, it must be set here.

See also:
[`dynamical_system`](@ref),
[`default_solver_options`](@ref),
[`default_parameters`](@ref),
[`solver_algorithm`](@ref),
"""
function FlowMap(sys;
                 solve_opts=DEFAULT_SOLVER_OPTIONS,
                 default_params=nothing,
                 algorithm=DEFAULT_INTEGRATION_ALGORITHM())
    FlowMap(
        sys,
        solve_opts,
        default_params,
        algorithm
    )
end

"""
    (::FlowMap)(q0, tspan; [p, solve_opts...])

Evaluate the mapping of `q0` in the flow map over time interval specified in `tspan`.

See also:
[`dynamical_system`](@ref),
[`default_solver_options`](@ref),
[`default_parameters`](@ref),
[`solver_algorithm`](@ref),
"""
function (fm::FlowMap)(q0, tspan; p::P=nothing, solve_opts...) where {P}
    if P == Nothing
        params = default_parameters(fm)
    else
        params = p
    end

    prob = ODEProblem(
        dynamical_system(fm),
        q0,
        tspan,
        params,
    )

    solve(
        prob,
        solver_algorithm(fm);
        default_solver_options(fm)...,
        solve_opts...
    )
end
