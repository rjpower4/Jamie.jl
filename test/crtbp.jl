@testset "CRTBP" begin
    @testset "CRTBP System" begin

        @testset "Parametric Typing" begin
            μ = 0.012
            name = "Example System"
            mstar = Float64(1.0)
            lstar = Float32(2.0)
            tstar = Float16(3.0)

            sys = CrtbpSystem(μ, name=name, char_mass=mstar, char_time=tstar, char_length=lstar)

            @test characteristic_mass(sys) == mstar
            @test typeof(characteristic_mass(sys)) == typeof(mstar)
            
            @test characteristic_time(sys) == tstar
            @test typeof(characteristic_time(sys)) == typeof(tstar)
            
            @test characteristic_length(sys) == lstar
            @test typeof(characteristic_length(sys)) == typeof(lstar)

            @test characteristic_velocity(sys) == (lstar / tstar)

            sys_nothing = CrtbpSystem(μ, name=name)
            @test characteristic_mass(sys_nothing) == nothing
            @test characteristic_length(sys_nothing) == nothing
            @test characteristic_time(sys_nothing) == nothing

        end
    end
end
