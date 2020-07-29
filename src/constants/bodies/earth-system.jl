struct Earth <: CelestialBody end
name(::Earth) = "EARTH"
spice_identifier(::Earth) = 399
equatorial_radius(::Earth) = 6378.137 # km
polar_radius(::Earth) = 6356.752 # km
mean_radius(::Earth) = 6371.000 # km
ellipticity(::Earth) = 0.003353
mean_density(::Earth) = 5514.0 # km / m^3
number_natural_satellites(::Earth) = 1
mean_semimajor_axis(::Earth) = 149.60e6 # km
mean_eccentricity(::Earth) = 0.0167086
solar_irradiance(::Earth) = 1361.0 # W / m^2
gravitational_parameter(::Earth) = 3.98600e5 # km^3 / s^2
total_mass(::Earth) = 5.9724e24 # kg
total_volume(::Earth) = 108.321e10 # km^3
parent_body(::Earth) = Sun()

struct Moon <: CelestialBody end

name(::Moon) = "MOON"
spice_identifier(::Moon) = 301
mean_radius(::Moon) = 1.7374e3 # km
gravitational_parameter(::Moon) = 4.9028e3 # km^3 / s^2
mean_semimajor_axis(::Moon) = 3.8474799201129237e5 # km
parent_body(::Moon) = Earth()
synodic_period(::Moon) = 29.53 # days