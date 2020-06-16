@testset "Bodies" begin
    earth_name = "EARTH"
    earth_spice_id = 399
    earth_potential = PointMassPotential(398600.4354360959)
    earth_shape = SphericalShapeModel(6371.008366666666)
    earth = CelestialBody(
        earth_name, 
        earth_spice_id, 
        earth_potential, 
        earth_shape
    ) 

    @test name_string(earth) == earth_name
    @test spice_identifier(earth) == earth_spice_id
    @test potential_model(earth) == earth_potential
    @test gravitational_parameter(earth) == gravitational_parameter(earth_potential)
    @test mean_radius(earth) == mean_radius(earth_shape)
    @test shape_model(earth) == earth_shape
end
