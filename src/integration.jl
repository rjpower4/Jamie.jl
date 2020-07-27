struct FlowMap{SYS, SO, P, A}
    system::SYS   # Dynamical System
    sopts::SO     # Default Solver Options
    params::P     # Default Parameters
    algorithm::A  # Solver algorithm for DifferentialEquations.jl
end

dynamical_system(fm::FlowMap) = fm.system

default_solver_options(fm::FlowMap) = fm.sopts

default_parameters(fm::FlowMap) = fm.params

solver_algorithm(fm::FlowMap) = fm.algorithm


function FlowMap(sys;
                 solve_opts=Constants.Integration.DEFAULT_SOLVE_OPTS,
                 default_params=nothing,
                 algorithm=Constants.Integration.DEFAULT_INTEGRATION_ALGORITHM())
    FlowMap(
        sys,
        solve_opts,
        default_params,
        algorithm
    )
end

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
