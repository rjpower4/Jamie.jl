@testset "Position-Velocity" begin
    x = 1.0
    y = 2.0
    z = 3.0
    vx = 4.0
    vy = 5.0
    vz = 6.0

    pos_mag = sqrt( x^2 +  y^2 +  z^2)
    vel_mag = sqrt(vx^2 + vy^2 + vz^2)


    pv1 = PositionVelocity(x, y, z, vx, vy, vz)

    @test pv1.x == x
    @test pv1.y == y
    @test pv1.z == z
    @test pv1.vx == vx
    @test pv1.vy == vy
    @test pv1.vz == vz

    @test pv_position(pv1) == SVector(x, y, z)
    @test pv_velocity(pv1) == SVector(vx, vy, vz)

    @test pv_position_mag(pv1) == pos_mag
    @test pv_velocity_mag(pv1) == vel_mag

    @test pv_position_unit(pv1) == (SVector(x, y, z) ./ pos_mag)
    @test pv_velocity_unit(pv1) == (SVector(vx, vy, vz) ./ vel_mag)
end
