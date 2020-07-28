const GRAVITATIONAL_CONSTANT = 6.6743e-20 # km^3 / (kg s^2)
const SPEED_OF_LIGHT = 299792.458 # km / s

include("integration.jl")
include("bodies/bodies.jl")

const EARTH_MOON_SYS = CrtbpSystem(
    1.2150584269940356e-02;
    name = "EARTH-MOON",
    dimset = DimensionalSet(
        mass = 6.0460429902763588e+24,
        length = 3.8474799201129237e+05,
        time = 3.7569985908499121e+05
    )
)

const LUNAR_SYNODIC_PERIOD = 29.530575 * (24.0 * 3600.0)