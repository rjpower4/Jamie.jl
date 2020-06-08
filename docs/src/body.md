# Celestial Bodies

## Potential Models
```@docs
NullPotential
PointMassPotential
gravitational_parameter(::PointMassPotential)
```

## Shape Models
```@docs
NullShapeModel
SphericalShapeModel
mean_radius(::SphericalShapeModel)
```

## The Celestial Body Structure
```@docs
CelestialBody
name_string(::CelestialBody)
spice_identifier
potential_model
shape_model
mean_radius(::CelestialBody)
gravitational_parameter(::CelestialBody)
```
