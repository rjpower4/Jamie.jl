module Constants

using Jamie: CelestialBody, PointMassPotential, SphericalShapeModel, CrtbpSystem

const EARTH = CelestialBody(
    "EARTH",
    399,
    PointMassPotential(3.9860043543609593e+05),
    SphericalShapeModel(6.3710083666666660e+03)
)

const MOON = CelestialBody(
    "MOON",
    301,
    PointMassPotential(4.9028000661637961e+03),
    SphericalShapeModel(1.7374000000000003e+03)
)

const EARTH_MOON_SYS = CrtbpSystem(
    1.2150584269940356e-02;
    name = "EARTH-MOON",
    char_mass = 6.0460429902763588e+24,
    char_length = 3.8474799201129237e+05,
    char_time = 3.7569985908499121e+05
)

const LUNAR_SYNODIC_PERIOD = 29.530575 * (24.0 * 3600.0)

end
