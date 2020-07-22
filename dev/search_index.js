var documenterSearchIndex = {"docs":
[{"location":"posvel.html#Position-Velocity-1","page":"Position-Velocity","title":"Position-Velocity","text":"","category":"section"},{"location":"posvel.html#","page":"Position-Velocity","title":"Position-Velocity","text":"PositionVelocity","category":"page"},{"location":"posvel.html#Jamie.PositionVelocity","page":"Position-Velocity","title":"Jamie.PositionVelocity","text":"PositionVelocity{T} <: FieldVector{6, T}\n\nStruct representing the position-velocity vector. T designates the types of component in the vector.\n\nSee also: pv_position, pv_velocity, pv_position_mag, pv_position_unit, pv_velocity_mag, pv_velocity_unit\n\n\n\n\n\n","category":"type"},{"location":"posvel.html#","page":"Position-Velocity","title":"Position-Velocity","text":"pv_position\npv_position_mag\npv_position_unit\npv_velocity\npv_velocity_mag\npv_velocity_unit","category":"page"},{"location":"posvel.html#Jamie.pv_position","page":"Position-Velocity","title":"Jamie.pv_position","text":"pv_position(::PositionVelocity)\n\nGet the position components of the position-velocity vector.\n\nSee also: pv_velocity, pv_position_mag, pv_position_unit, pv_velocity_mag, pv_velocity_unit\n\n\n\n\n\n","category":"function"},{"location":"posvel.html#Jamie.pv_position_mag","page":"Position-Velocity","title":"Jamie.pv_position_mag","text":"pv_position_mag(::PositionVelocity)\n\nReturn the magnitude of the position components of the position-velocity vector. This is equivalent to the distance from the origin of the position vector.\n\nSee also: pv_position, pv_velocity, pv_position_unit, pv_velocity_mag, pv_velocity_unit\n\n\n\n\n\n","category":"function"},{"location":"posvel.html#Jamie.pv_position_unit","page":"Position-Velocity","title":"Jamie.pv_position_unit","text":"pv_position_unit(::PositionVelocity)\n\nReturn the unit vector parallel to the position sub-vector in the specified position-velocity vector.\n\nSee also: pv_position, pv_velocity, pv_position_mag, pv_velocity_mag, pv_velocity_unit\n\n\n\n\n\n","category":"function"},{"location":"posvel.html#Jamie.pv_velocity","page":"Position-Velocity","title":"Jamie.pv_velocity","text":"pv_velocity(::PositionVelocity)\n\nGet the velocity components of the position-velocity vector.\n\nSee also: pv_position, pv_position_mag, pv_position_unit, pv_velocity_mag, pv_velocity_unit\n\n\n\n\n\n","category":"function"},{"location":"posvel.html#Jamie.pv_velocity_mag","page":"Position-Velocity","title":"Jamie.pv_velocity_mag","text":"pv_velocity_mag(::PositionVelocity)\n\nReturn the magnitude of the velocity components of the position-velocity vector.\n\nSee also: pv_position, pv_velocity, pv_position_mag, pv_position_unit, pv_velocity_unit\n\n\n\n\n\n","category":"function"},{"location":"posvel.html#Jamie.pv_velocity_unit","page":"Position-Velocity","title":"Jamie.pv_velocity_unit","text":"pv_velocity_unit(::PositionVelocity)\n\nReturn the unit vector parallel to the velocity sub-vector in the specified position-velocity vector.\n\nSee also: pv_position, pv_velocity, pv_position_mag, pv_position_unit, pv_velocity_mag,\n\n\n\n\n\n","category":"function"},{"location":"body.html#Celestial-Bodies-1","page":"Celestial Bodies","title":"Celestial Bodies","text":"","category":"section"},{"location":"body.html#Potential-Models-1","page":"Celestial Bodies","title":"Potential Models","text":"","category":"section"},{"location":"body.html#","page":"Celestial Bodies","title":"Celestial Bodies","text":"NullPotential\nPointMassPotential\ngravitational_parameter(::PointMassPotential)","category":"page"},{"location":"body.html#Jamie.NullPotential","page":"Celestial Bodies","title":"Jamie.NullPotential","text":"NullPotential\n\nEmpty potential model used as an indicator that no potential is defined.\n\nSee also: PointMassPotential\n\n\n\n\n\n","category":"type"},{"location":"body.html#Jamie.PointMassPotential","page":"Celestial Bodies","title":"Jamie.PointMassPotential","text":"PointMassPotential(gm::T)\n\nSimple potential model for a point mass (or centrobaric) body of gravitational parameter gm.\n\nSee also: NullPotential\n\n\n\n\n\n","category":"type"},{"location":"body.html#Jamie.gravitational_parameter-Tuple{PointMassPotential}","page":"Celestial Bodies","title":"Jamie.gravitational_parameter","text":"gravitational_parameter(::PointMassPotential)\n\nRetrieve the gravitational parameter, GM, of the potential model. This is the standard 2-body μ value.\n\nSee also: PointMassPotential\n\n\n\n\n\n","category":"method"},{"location":"body.html#Shape-Models-1","page":"Celestial Bodies","title":"Shape Models","text":"","category":"section"},{"location":"body.html#","page":"Celestial Bodies","title":"Celestial Bodies","text":"NullShapeModel\nSphericalShapeModel\nmean_radius(::SphericalShapeModel)","category":"page"},{"location":"body.html#Jamie.NullShapeModel","page":"Celestial Bodies","title":"Jamie.NullShapeModel","text":"NullShapeModel\n\nEmpty shape model used as an indicator that no shape is defined.\n\nSee also: SphericalShapeModel\n\n\n\n\n\n","category":"type"},{"location":"body.html#Jamie.SphericalShapeModel","page":"Celestial Bodies","title":"Jamie.SphericalShapeModel","text":"SphericalShapeModel(radius::T)\n\nSimple shape model describing a sphere of radius, radius.\n\nSee also: NullShapeModel\n\n\n\n\n\n","category":"type"},{"location":"body.html#Jamie.mean_radius-Tuple{SphericalShapeModel}","page":"Celestial Bodies","title":"Jamie.mean_radius","text":"mean_radius(::SphericalShapeModel)\n\nRetrieve the mean radius of the spherical shape model. Clearly for this model, the mean radius is simple equal to the radius, radius.\n\nSee also: NullShapeModel\n\n\n\n\n\n","category":"method"},{"location":"body.html#The-Celestial-Body-Structure-1","page":"Celestial Bodies","title":"The Celestial Body Structure","text":"","category":"section"},{"location":"body.html#","page":"Celestial Bodies","title":"Celestial Bodies","text":"CelestialBody\nname_string(::CelestialBody)\nspice_identifier\npotential_model\nshape_model\nmean_radius(::CelestialBody)\ngravitational_parameter(::CelestialBody)","category":"page"},{"location":"body.html#Jamie.CelestialBody","page":"Celestial Bodies","title":"Jamie.CelestialBody","text":"CelestialBody{T, I <: Integer}(name::String, gm::T, spice_id::I, mean_radius::T) <: AbstractCelestialBody\n\nStruct defining a celestial body (e.g. Moon, Planet, Asteroid).\n\nSee also: PointMassPotential, SphericalShapeModel, NullCelestialBody\n\n\n\n\n\n","category":"type"},{"location":"body.html#Jamie.name_string-Tuple{CelestialBody}","page":"Celestial Bodies","title":"Jamie.name_string","text":"name_string(cb::CelestialBody)\n\nRetrieve the name of the celestial body.\n\nSee also: CelestialBody, spice_identifier, potential_model, shape_model, mean_radius, gravitational_parameter\n\n\n\n\n\n","category":"method"},{"location":"body.html#Jamie.spice_identifier","page":"Celestial Bodies","title":"Jamie.spice_identifier","text":"spice_identifier(cb::CelestialBody)\n\nRetrieve the integer SPICE identifier of the celestial body. More information on SPICE ID codes can be found here.\n\nSee also: CelestialBody, name_string, potential_model, shape_model, mean_radius, gravitational_parameter\n\n\n\n\n\n","category":"function"},{"location":"body.html#Jamie.potential_model","page":"Celestial Bodies","title":"Jamie.potential_model","text":"potential_model(::NullCelestialBody)\n\nReturns a null potential model.\n\nSee also: NullPotentialModel\n\n\n\n\n\npotential_model(::CelestialBody)\n\nRetrieve the potential modlel of the body.\n\nSee also: CelestialBody, name_string, spice_identifier, shape_model, mean_radius, gravitational_parameter\n\n\n\n\n\n","category":"function"},{"location":"body.html#Jamie.shape_model","page":"Celestial Bodies","title":"Jamie.shape_model","text":"shape_model(::NullCelestialBody)\n\nReturns a null shape model.\n\nSee also: NullShapeModel\n\n\n\n\n\nshape_model(::CelestialBody)\n\nRetrieve the shape model of the body.\n\nSee also: CelestialBody, name_string, spice_identifier, potential_model, mean_radius, gravitational_parameter\n\n\n\n\n\n","category":"function"},{"location":"body.html#Jamie.mean_radius-Tuple{CelestialBody}","page":"Celestial Bodies","title":"Jamie.mean_radius","text":"mean_radius(::CelestialBody)\n\nRetrieve the mean radius of the body if a non-null shape model is defined. The particular calculation of the mean radius is implemented by the shape model.\n\nSee also: SphericalShapeModel, CelestialBody, name_string, spice_identifier, potential_model, shape_model, gravitational_parameter\n\n\n\n\n\n","category":"method"},{"location":"body.html#Jamie.gravitational_parameter-Tuple{CelestialBody}","page":"Celestial Bodies","title":"Jamie.gravitational_parameter","text":"gravitational_parameter(::CelestialBody)\n\nRetrieve the gravitational parameter of the body if a non-null potential model is defined. The particular calculation of the gravitational parameter is implemented by the potential model.\n\nSee also: PointMassPotential, CelestialBody, name_string, spice_identifier, potential_model, shape_model, mean_radius,\n\n\n\n\n\n","category":"method"},{"location":"crtbp.html#Circular-Restricted-Three-Body-Problem-(CRTBP)-1","page":"CRTBP","title":"Circular Restricted Three Body Problem (CRTBP)","text":"","category":"section"},{"location":"crtbp.html#Structures-1","page":"CRTBP","title":"Structures","text":"","category":"section"},{"location":"crtbp.html#","page":"CRTBP","title":"CRTBP","text":"CrtbpP1\nCrtbpP2\nCrtbpSystem","category":"page"},{"location":"crtbp.html#Jamie.CrtbpP1","page":"CRTBP","title":"Jamie.CrtbpP1","text":"CrtbpP1\n\nEmpty struct identifying the first (usually larger) primary in the CRTBP.\n\nSee also: CrtbpP2\n\n\n\n\n\n","category":"type"},{"location":"crtbp.html#Jamie.CrtbpP2","page":"CRTBP","title":"Jamie.CrtbpP2","text":"CrtbpP2\n\nEmpty struct identifying the second (usually smaller) primary in the CRTBP.\n\nSee also: CrtbpP1\n\n\n\n\n\n","category":"type"},{"location":"crtbp.html#Jamie.CrtbpSystem","page":"CRTBP","title":"Jamie.CrtbpSystem","text":"CrtbpSystem{T, CM, CL, CT}(μ::T, name::String, char_mass::CM, char_length::CL, char_time::CT)\n\nStruct defining the parameters for the Crtbp System. Contains the mass ratio parameter, mu, as well as the system name and characteristic quantities.\n\nSee also: mass_ratio, characteristic_mass, characteristic_time, characteristic_length, characteristic_velocity\n\n\n\n\n\n","category":"type"},{"location":"crtbp.html#Functions-1","page":"CRTBP","title":"Functions","text":"","category":"section"},{"location":"crtbp.html#","page":"CRTBP","title":"CRTBP","text":"mass_ratio\ncharacteristic_mass\ncharacteristic_length\ncharacteristic_time\ncharacteristic_velocity\ndimensionalize\nnondimensionalize\npseudopotential\npseudopotential_gradient\npseudopotential_hessian\njacobi_constant\njacobi_constant_gradient\ncrtbp_eom\ncrtbp_jacobian","category":"page"},{"location":"crtbp.html#Jamie.mass_ratio","page":"CRTBP","title":"Jamie.mass_ratio","text":"mass_ratio(::CrtbpSystem)\n\nRetrieve the mass ratio, μ, of the specified system.\n\nSee also: CrtbpSystem, characteristic_mass, characteristic_time, characteristic_length, characteristic_velocity\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.characteristic_mass","page":"CRTBP","title":"Jamie.characteristic_mass","text":"characteristic_mass(::DimensionalSet)\n\nRetrieve the characteristic mass of the dimensional set specified. \n\nSee also:  DimensionalSet, characteristic_time, characteristic_length,  characteristic_velocity, characteristic_acceleration\n\n\n\n\n\ncharacteristic_mass(::CrtbpSystem)\n\nRetrieve the characteristic mass of the CRTBP system specified. This mass is usually defined as the sum of the masses of the two primaries.\n\nSee also: CrtbpSystem,  mass_ratio,  dimensional_set, characteristic_time, characteristic_length,  characteristic_velocity, characteristic_acceleration\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.characteristic_length","page":"CRTBP","title":"Jamie.characteristic_length","text":"characteristic_length(::DimensionalSet)\n\nRetrieve the characteristic length of the dimensional set specified. \n\nSee also:  DimensionalSet, characteristic_mass,  characteristic_time, characteristic_velocity, characteristic_acceleration\n\n\n\n\n\ncharacteristic_length(::CrtbpSystem)\n\nRetrieve the characteristic length of the CRTBP system specified. This length is usually defined as the semi-major axis of the circular orbit of P2 in the CRTBP.\n\nSee also: CrtbpSystem,  dimensional_set, mass_ratio, characteristic_mass, characteristic_time,  characteristic_acceleration\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.characteristic_time","page":"CRTBP","title":"Jamie.characteristic_time","text":"characteristic_time(::DimensionalSet)\n\nRetrieve the characteristic time of the dimensional set specified. \n\nSee also:  DimensionalSet, characteristic_mass,  characteristic_length, characteristic_velocity, characteristic_acceleration\n\n\n\n\n\ncharacteristic_length(::CrtbpSystem)\n\nRetrieve the characteristic time of the CRTBP system specified. This time is usually defined as the time such that the universal gravitational constant in the non-dimensional system is unity.\n\nSee also: CrtbpSystem,  dimensional_set, mass_ratio, characteristic_mass, characteristic_length,  characteristic_velocity, characteristic_acceleration\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.characteristic_velocity","page":"CRTBP","title":"Jamie.characteristic_velocity","text":"characteristic_velocity(::DimensionalSet)\n\nRetrieve the characteristic velocity of the dimensional set specified.  Note, this is a derived quantity and equal to the characteristic length over the  characteristic time.\n\nSee also:  DimensionalSet, characteristic_mass,  characteristic_length, characteristic_time, characteristic_acceleration\n\n\n\n\n\ncharacteristic_velocity(::CrtbpSystem)\n\nRetrieve the characteristic velocityof the CRTBP system specified. This is equivalent to the characteristic length divided by the characteristic time.\n\nSee also: CrtbpSystem,  dimensional_set, mass_ratio, characteristic_mass, characteristic_length, characteristic_acceleration\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.dimensionalize","page":"CRTBP","title":"Jamie.dimensionalize","text":"dimensionalize(::CrtbpSystem{T,T}, ::PositionVelocity)\n\nDimensionalize the position and velocity components using the characteristic quantities associated with the CRTBP system.\n\nSee also: nondimensionalize CrtbpSystem, characteristic_length, characteristic_velocity\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.nondimensionalize","page":"CRTBP","title":"Jamie.nondimensionalize","text":"nondimensionalize(::CrtbpSystem{T,T}, ::PositionVelocity)\n\nNondimensionalize the position and velocity components using the characteristic quantities associated with the CRTBP system.\n\nSee also: dimensionalize CrtbpSystem, characteristic_length, characteristic_velocity\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.pseudopotential","page":"CRTBP","title":"Jamie.pseudopotential","text":"pseudopotential(::CrtbpSystem, ::PositionVelocity)\n\nCalculate the CRTBP pseudopotential of the specified position-velocity vector.\n\nSee also: pseudopotential_gradient, pseudopotential_hessian\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.pseudopotential_gradient","page":"CRTBP","title":"Jamie.pseudopotential_gradient","text":"pseudopotential_gradient(::CrtbpSystem, ::PositionVelocity)\n\nCalculate the partials of the CRTBP pseudopotential with respect to all 6 state components in the position-velocity vector.\n\nSee also: pseudopotential, pseudopotential_hessian\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.pseudopotential_hessian","page":"CRTBP","title":"Jamie.pseudopotential_hessian","text":"pseudopotential_hessian(::CrtbpSystem, ::PositionVelocity)\n\nCalculate the second order partials of the CRTBP pseudopotential with respect to all 6 state components in the position-velocity vector.\n\nSee also: pseudopotential, pseudopotential_gradient\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.jacobi_constant","page":"CRTBP","title":"Jamie.jacobi_constant","text":"jacobi_constant(::CrtbpSystem, ::PositionVelocity)\n\nCalculate the Jacobi constant of a state in the CRTBP.\n\nSee also: jacobi_constant_gradient\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.jacobi_constant_gradient","page":"CRTBP","title":"Jamie.jacobi_constant_gradient","text":"jacobi_constant_gradient(::CrtbpSystem, ::PositionVelocity)\n\nCalculate the gradient of Jacobi constant with respect to all 6 position-velocity variables.\n\nSee also: jacobi_constant\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.crtbp_eom","page":"CRTBP","title":"Jamie.crtbp_eom","text":"crtbp_eom(pv::AbstractArray, sys::CrtbpSystem, t)\n\nEvaluate the equations of motion for the Circular Restricted Three Body Problem (CRTBP) at the state specified by pv in the system specified by sys. The time input, t, while required due to compatibility with DifferentialEquations.jl, is not used.\n\nNote, here, pv can be any AbstractArray due to compatibility with DifferentialEquations. However, if a PositionVelocity is not constructible via PositionVelocity(pv), an error will be thrown.\n\nSee also: crtbp_jacobian, PositionVelocity\n\n\n\n\n\n","category":"function"},{"location":"crtbp.html#Jamie.crtbp_jacobian","page":"CRTBP","title":"Jamie.crtbp_jacobian","text":"crtbp_jacobian(pv::AbstractArray, sys::CrtbpSystem, t)\n\nEvaluate the jacobian of equations of motion for the Circular Restricted Three Body Problem (CRTBP) at the state specified by pv in the system specified by sys. The time input, t, while required due to compatibility with DifferentialEquations.jl, is not used.\n\nNote, here, pv can be any AbstractArray due to compatibility with DifferentialEquations. However, if a PositionVelocity is not constructible via PositionVelocity(pv), an error will be thrown.\n\nSee also: crtbp_eom, PositionVelocity\n\n\n\n\n\n","category":"function"},{"location":"index.html#Jamie.jl-Documentation-1","page":"Introduction","title":"Jamie.jl Documentation","text":"","category":"section"},{"location":"util.html#Utilities-1","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"util.html#Generating-Points-on-a-Sphere-1","page":"Utilities","title":"Generating Points on a Sphere","text":"","category":"section"},{"location":"util.html#","page":"Utilities","title":"Utilities","text":"points_on_sphere","category":"page"},{"location":"util.html#Jamie.points_on_sphere","page":"Utilities","title":"Jamie.points_on_sphere","text":"points_on_sphere(n; [radius=1.0], [center=(0.0, 0.0, 0.0)])\n\nGenerate a set of points on a sphere approximately uniformly distributed via the Fibonacci lattice method. The points are returned as an Array of n SVectors of length 3. A specific radius and center can be optionally given as key-word arguments.\n\n\n\n\n\n","category":"function"}]
}
