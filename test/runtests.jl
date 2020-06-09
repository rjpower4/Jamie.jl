using Jamie
using Test
using StaticArrays

# Hook into Pkg.test so that tests from a single file can be run.  For example,
# to run only the MVector and SVector tests, use:
#
#   Pkg.test("StaticArrays", test_args=["MVector", "SVector"])
#
# Source: https://github.com/JuliaArrays/StaticArrays.jl/blob/master/test/runtests.jl
enabled_tests = lowercase.(ARGS)
function addtests(fname)
    key = lowercase(splitext(fname)[1])
    if isempty(enabled_tests) || key in enabled_tests
        include(fname)
    end
end


addtests("posvel.jl")
addtests("crtbp.jl")
addtests("body.jl")
