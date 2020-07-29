@testset "Bodies" begin
    earth_name = "EARTH"
    earth_spice_id = 399
    earth = Jamie.Earth()

    @test name(earth) == earth_name
    @test spice_identifier(earth) == earth_spice_id
end
